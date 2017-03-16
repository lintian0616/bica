import os, sys
import string

def XenomeCorrect(XenomeOutputFastq):
    Corrected = ""
    with open(XenomeOutputFastq, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count % 4 == 1:
                line = "@" + line
            elif count % 4 == 3:
                line = "+" + line
            Corrected = Corrected + line
    Corrected = Corrected[:-1]
    print Corrected

if __name__ == '__main__':
    filename = sys.argv[1]
    XenomeCorrect(filename)
