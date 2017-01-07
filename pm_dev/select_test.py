import sys
import select

while(True):
    i,o,e = select.select([sys.stdin.read(1)],[],[],10)
    print(i)
