import ConfigParser
import sys
import os


def createDataPath(dataPath='data'):
    """Returns the Current Path of the "executable" """
    #Check if running from the binary
    frozen = getattr(sys, 'frozen', None)
    if frozen:
        programPath= os.path.split(__file__.replace("%slibrary.zip"%os.path.sep,""))[0]
    else:
        programPath= os.path.split(__file__)[0]
    return os.path.join(programPath, dataPath)


def ReadConfig(configFile="SupportConfig.cfg"):
    config=ConfigParser.ConfigParser()
    config.read(configFile)
    return config

if __name__ =="__main__":
    config=ReadConfig()
    
    print config.sections()
    print config.getfloat('parameters','window')
    print config.getfloat('parameters','generations')

    print config.items('data')
    
    
    
