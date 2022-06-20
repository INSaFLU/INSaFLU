'''
Created on 20/06/2022

@author: mmp
'''

class CmdException(Exception):
    '''
    classdocs
    '''


    def __init__(self, message, cmd, output_path = ""):
        '''
        Constructor
        '''
        self.message = message
        self.cmd = cmd
        self.output_path = output_path
        super().__init__(self.message)
        
    def __str__(self):
        return "{}\n{}".format(self.message, self.cmd)
        
    def exist_path(self):
        return len(self.output_path) > 0