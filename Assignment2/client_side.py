import multiprocessing as mp
from multiprocessing.managers import BaseManager
import queue
import time


__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "14-5-2022"


class ClientSide:
    """
    Class that handles the client side
    """

    def __init__(self, ip: str, port: int, auth_key: str, poison_pill: str) -> None:
        self.ip = ip
        self.port = port
        self.auth_key = auth_key
        self.poison_pill = poison_pill
        self.error = "DOH"

    def make_client_manager(self) -> BaseManager:
        """ 
        Create a manager for a client. This manager connects to a server on the
        given address and exposes the get_job_q and get_result_q methods for 
        accessing the shared queues from the server.
            
        :returns
        --------
        manager - ServerQueueManager
            manager object
        """
        class ServerQueueManager(BaseManager):
            pass

        ServerQueueManager.register('get_job_q')
        ServerQueueManager.register('get_result_q')

        manager = ServerQueueManager(address=(self.ip, self.port), authkey=self.auth_key)
        manager.connect()

        print('Client connected to %s - %s' % (self.ip, self.port))
        return manager

    def run_client(self, num_processes: int) -> None:
        """
        Run a client by making a client manager and starting the workers.

        :parameters
        -----------
        num_processes - int
            Number of processes to start (number of peons to start)
        """
        manager = self.make_client_manager()
        job_q = manager.get_job_q()
        result_q = manager.get_result_q()
        self.run_workers(job_q, result_q, num_processes)
        
    def run_workers(self, job_q: queue, result_q: queue, num_processes: int) -> None:
        """
        Run the workers by starting n peons. 

        :parameters
        -----------
        job_q - queue
            Queue containing the jobs to do
        result_q - queue
            Queue containing the results
        num_processes - int
            Number of processes to start (number of peons to start) 
        """
        processes = []
        for p in range(num_processes):
            temP = mp.Process(target=self.peon, args=(job_q, result_q))
            processes.append(temP)
            temP.start()
        print("Started %s workers!" % len(processes))
        for temP in processes:
            temP.join()

    def peon(self, job_q: queue, result_q: queue) -> None:
        """
        Process the data and store the result in the result queue.
        Keep processing data untill the poisonpill is found.

        :parameters
        -----------
        job_q - queue
            Queue containing the jobs to do
        result_q - queue
            Queue containing the results
        """
        my_name = mp.current_process().name
        while True:
            try:
                job = job_q.get_nowait()
                if job == self.poison_pill:
                    job_q.put(self.poison_pill)
                    print("Aaaaaaargh", my_name)
                    return
                else:
                    try:
                        result = job['fn'](job['arg'])
                        print("Peon %s Workwork on %s!" % (my_name, job['arg']))
                        result_q.put({'job': job, 'result' : result})
                    except NameError:
                        print("Can't find yer fun Bob!")
                        result_q.put({'job': job, 'result' : self.error})

            except queue.Empty:
                print("sleepytime for", my_name)
                time.sleep(1)