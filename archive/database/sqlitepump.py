#!/usr/bin/env python

import os, sys, tempfile, csv
import re, inspect, ast
from subprocess import Popen, PIPE, STDOUT

"""
Class for dumping data into a SQLite database.
You will need the sqlite3 command line client in your PATH,
even if you are on Windows and not UNIX (maybe)

The simplest use case is as follows:

# import the class
from sqlitemagic import SQLiteMagician

# Your data should be a list of equal length tuples
data = [(x,) for x in range(100)]

# Create an SQLiteMagician to do stuff for you
wizzard = SQLiteMagician('my_database.db')

# Conjure a table with all your data populated 
# and the table created with the same name as 
# the variable passed in, so 'data' in this case.
wizzard.conjure(data)

2013-08-10 Matt Oates
MattOates[AT]gmail[DOT]com
>:3
"""
class SQLiteMagician(object):

    #Default command line client to use to access SQLite3
    command = "sqlite3"

    #Python to SQLite type map, anything not in here will be mapped to NONE
    types = {"int":"INTEGER","str":"TEXT","float":"REAL"}

    #Constructor, what db to use 
    def __init__(self,database):
        self.database = database
    
    #Using an example data record infer the SQLite field types.
    #Field names will be made up if a header list was not specified
    def _get_field_spec(self,data,header=None,fieldname="field"):
        types = [
                "NONE" if t == "NoneType" else t 
                    for t in [ SQLiteMagician.types.get( type( d ).__name__ ) 
                        for d in data[0] 
                    ] 
                ]

        #If no header was passed in create a load of dumby field titles
        if header == None:
            header = ["field%s" % n for n in range(len(data[0]))]
        
        #Return the fieldspec e.g. "field1 INTEGER, field2 TEXT"
        return ", ".join(["%s %s" % field for field in zip(header,types)])
        
    #Convoluted frame stack walk and source scrape to get what the calling statement to a function looked like.
    #Specifically return the name of the variable passed as parameter found at position pos in the parameter list.
    def _caller_param_name(self, pos):
        #The parameter name to return
        param = None
        #Get the frame object for this function call
        thisframe = inspect.currentframe()
        try:
            #Get the parent calling frames details
            frames = inspect.getouterframes(thisframe)
            #Function this function was just called from that we wish to find the calling parameter name for
            function = frames[1][3]
            #Get all the details of where the calling statement was
            frame,filename,line_number,function_name,source,source_index = frames[2]
            #Read in the source file in the parent calling frame upto where the call was made
            with open(filename) as source_file:
                head=[source_file.next() for x in xrange(line_number)]
            source_file.close()
    
            #Build all lines of the calling statement, this deals with when a function is called with parameters listed on each line
            lines = []
            #Compile a regex for matching the start of the function being called
            regex = re.compile(r'\.?\s*%s\s*\(' % (function))
            #Work backwards from the parent calling frame line number until we see the start of the calling statement (usually the same line!!!)
            for line in reversed(head):
                lines.append(line.strip())
                if re.search(regex, line):
                    break
            #Put the lines we have groked back into sourcefile order rather than reverse order
            lines.reverse()
            #Join all the lines that were part of the calling statement
            call = "".join(lines)
            #Grab the parameter list from the calling statement for the function we were called from
            match = re.search('\.?\s*%s\s*\((.*)\)' % (function), call)
            paramlist = match.group(1)
            #If the function was called with no parameters raise an exception
            if paramlist == "":
                raise LookupError("Function called with no parameters.")
            #Use the Python abstract syntax tree parser to create a parsed form of the function parameter list 'Name' nodes are variable names
            parameter = ast.parse(paramlist).body[0].value
            #If there were multiple parameters get the positional requested
            if type(parameter).__name__ == 'Tuple':
                #If we asked for a parameter outside of what was passed complain
                if pos >= len(parameter.elts):
                    raise LookupError("The function call did not have a parameter at postion %s" % pos)
                parameter = parameter.elts[pos]
            #If there was only a single parameter and another was requested raise an exception
            elif pos != 0:
                raise LookupError("There was only a single calling parameter found. Parameter indices start at 0.")
            #If the parameter was the name of a variable we can use it otherwise pass back None
            if type(parameter).__name__ == 'Name':
                param = parameter.id
        finally:
            #Remove the frame reference to prevent cyclic references screwing the garbage collector
            del thisframe
        #Return the parameter name we found
        return param

    #Create a table given a single record of data and perhaps a header of field names
    def create_table(self,table,data,header=None,header_included=False):

        #Grab the header if its the first record of the data
        if header_included:
            header = data.pop(0)
        
        #Get the field spec for this table
        field_spec = self._get_field_spec(data,header)

        #Run the query using a subprocess to the sqlite3 command
        sqlite3 = Popen([SQLiteMagician.command, self.database], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        db = sqlite3.communicate(input='CREATE TABLE IF NOT EXISTS %s(%s);\n.quit\n' % (table,field_spec))[0]
    
    #Dump a variable to an SQLite database using the maximum amount of magic possible
    def conjure(self,data,table=None,header=None,header_included=False):
        #If the table name was not specified use stack inspection to see what the data variable was called!
        #Magical and horrible, but I want an interface thats incredibly simple for the simplest case of dumping a variable.
        if not table:
            table = self._caller_param_name(0)
        #Create a table to put the data into
        self.create_table(table,data,header,header_included)
        #Dump as normal
        self.dump(data,table)

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
            sqlite3 = Popen([SQLiteMagician.command, self.database], stdout=PIPE, stdin=PIPE, stderr=STDOUT)

            #Tell sqlite3 to import from our named pipe
            db = sqlite3.communicate(input='.separator "\\t"\n.import %s %s\n.quit\n' % (data_pipe,table))[0]
            sys.exit(1)

        #If we are the parent process stream data to the child
        else:
            #File handle to pipe to write data into table
            data_pipe_write_handle = open(data_pipe,'w')
            #Have the original program send data down the pipe
            writer = csv.writer(data_pipe_write_handle, delimiter="\t")
            writer.writerows(data)
            #Close the pipe once we are done!
            data_pipe_write_handle.close()
            #Wait for SQLite3 to finish importing, waiting on all child processes
            os.wait()
            #Remove the named pipe file we created because its junk and we dont want a clash
            os.unlink(data_pipe)
