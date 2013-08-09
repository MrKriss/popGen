'''
Created on 9 Aug 2013

@author: musselle
'''

#!/usr/bin/env python
 
import os, sys, tempfile, csv
from subprocess import Popen, PIPE, STDOUT
 

class SQLitePumper(object):
 
    command = 'sqlite3'
 
    def __init__(self,database):
        self.database = database
 
    def dump(self,data,table):
 
        #Name of the pipe, use tempfile to create some random filename usually in /tmp
        data_pipe = tempfile.mktemp('datapump')
        #Create the actual pipe 'file' where the OS knows you wanted a pipe type thing
        os.mkfifo( data_pipe, 0644 )
 
        #Create a child process to run in parallel
        child = os.fork()
 
        #If child is 0 we are the child
        if child == 0:
            #Have the child send the command to sqlite3
            #Create a forked process running sqlite3 with a pipe that we can send commands to
            sqlite3 = Popen([self.command, self.database], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
 
            #Tell sqlite3 to import from our named pipe
            db = sqlite3.communicate(input='.separator "\\t"\n.import %s %s\n.quit\n' % (data_pipe,table))[0]
            
            #The child exits so we stop waiting
            sys.exit(1)
 
        else:
            #File handle to pipe to write data into table
            data_pipe_write_handle = open(data_pipe,'w')
            #Send data down the pipe
            writer = csv.writer(data_pipe_write_handle, delimiter="\t")
            writer.writerows(data)
            data_pipe_write_handle.close()
            #Wait for SQLite3 to finish importing, waiting on all child processes
            os.wait()
            #Remove the named pipe file we created because its junk and we dont want a clash
            os.unlink(data_pipe)
            
            
if __name__ == '__main__':
    
#!/usr/bin/env python

from datapump import SQLitePumper

    pump = SQLitePumper( 'test.db' )
    data = [ (x,) for x in range(100) ]
    pump.dump( data, 'test' )
            