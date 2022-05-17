from multiprocessing.managers import BaseManager, SyncManager
import queue

class ProxyManager:

    def __init__(self, ip: str, port: int, auth_key: str):
        self.ip = ip
        self.port = port
        self.auth_key = auth_key

    def make_server_manager(self):
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
            pass

        QueueManager.register('get_job_q', callable=lambda: job_q)
        QueueManager.register('get_result_q', callable=lambda: result_q)

        manager = QueueManager(address=(self.ip, self.port), authkey=self.auth_key)
        manager.start()
        print('Server started at port %s' % self.port)
        return manager


    def make_client_manager(self):
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

        print('Client connected to %s:%s' % (self.ip, self.port))
        return manager