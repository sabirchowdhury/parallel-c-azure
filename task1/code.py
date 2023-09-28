import time
import os

numberOfThreads = [thread for thread in range(51,0,-10)]

with open("main.csv", "a") as csvFile:
    csvFile.write("\nArray size: 10000\nPrecision: 0.01\n\nThread,")
    for thread in numberOfThreads:
        csvFile.write(str(thread)+",")

    csvFile.write("\nTime(s),")
    csvFile.close()

    for thread in numberOfThreads:
        os.system("fahim -size 10000 -precision 0.01 -nthreads -border 1 -fill 0 %d >> main.csv"%thread)



