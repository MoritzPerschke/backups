import subprocess
import statistics as stats
import csv
runPrograms = []
dataStructure = [0, 1]
instructionMix = [0, 1, 10, 50]
elementSize = 8
numberOfElements = [10, 10000000]

for ds in dataStructure:
    for insMix in instructionMix:
        for numbElem in numberOfElements:
            runPrograms.append("./benchmark " + str(ds) + " " + str(insMix) + " " + str(elementSize) + " " + str(numbElem))

dictResults = {}
print(runPrograms)

#initialize dict
for k in runPrograms:
    dictResults[k] = []
for i in range (0,5):
    for k in runPrograms:
        result = ['Run' + str(i)]
        command = k
        proc = subprocess.run([command], shell=True, stdout=subprocess.PIPE)
        lastElement = proc.stdout.splitlines()[-1]
        dictResults[k].append(int(str(lastElement, 'UTF-8')))
       

for k in runPrograms:
    sum = 0
    for j in dictResults[k]:
        sum = sum + j
    sum = sum / 5
    sum = [sum]
    print(sum)
    with open(k + '.csv', 'w', newline='') as file: 
        writer = csv.writer(file)
        writer.writerow(sum)