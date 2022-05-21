"""
This module contains code from the server side.
"""

from multiprocessing.managers import BaseManager
import queue
import time


__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "14-5-2022"


class ServerSide:
    """
    Class that handles the server side.
    """

    def __init__(self, ip_adress: str, port: int, auth_key: str, poison_pill: str) -> None:
        self.ip_adress = ip_adress
        self.port = port
        self.auth_key = auth_key
        self.poison_pill = poison_pill

    def make_server_manager(self) -> BaseManager:
        """
        Create a manager for the server, listening on the given port.

        :returns
        --------
        manager - QueueManager
            manager object with get_job_q and get_result_q methods.
        """
        job_q = queue.Queue()
        result_q = queue.Queue()

        class QueueManager(BaseManager):
            """
            Server Queue Manager
            """

        QueueManager.register('get_job_q', callable=lambda: job_q)
        QueueManager.register('get_result_q', callable=lambda: result_q)

        manager = QueueManager(address=(self.ip_adress, self.port), authkey=self.auth_key)
        manager.start()
        print(f'Server started at port {self.port}')
        return manager


    def run_server(self, func_name, data) -> None:
        """
        Put jobs in the jobs queue and check if the data is being processed
        by checking if the result queue is getting filled. Stop all the clients when
        all data has been processed.

        :parameters
        -----------
        func_name - function
            Name of the function that will process the data
        data - array like
            Data that needs to be processed
        """
        # Start a shared manager server and access its queues
        manager = self.make_server_manager()
        shared_job_q = manager.get_job_q()
        shared_result_q = manager.get_result_q()

        if not data:
            print("Gimme something to do here!")
            return

        print("Sending data!")
        for dat in data:
            shared_job_q.put({'func_name' : func_name, 'arg' : dat})

        time.sleep(2)

        results = []
        while True:
            try:
                result = shared_result_q.get_nowait()
                results.append(result)
                print("Got result!", result)
                if len(results) == len(data):
                    print("Got all results!")
                    break
            except queue.Empty:
                time.sleep(1)
                continue
        # Tell the client process no more data will be forthcoming
        print("Time to kill some peons!")
        shared_job_q.put(self.poison_pill)
        # Sleep a bit before shutting down the server - to give clients time to
        # realize the job queue is empty and exit in an orderly way.
        time.sleep(5)
        print("Aaaaaand we're done for the server!")
        manager.shutdown()
        print(results)
