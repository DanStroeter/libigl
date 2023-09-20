import sys
import os
import subprocess
import csv

coordinateFlags = ['--MVC', '--harmonic', '--green', '--MEC', '--BBW', '--LBC', '--QGC', '--MLC', '--somigliana']

# In matching order
models = ['chessBishop.obj']
cages = ['bishop_cages_triangulated.off']
cagesDeformed = ['bishop_cages_triangulated_deformed.obj']
embeddings = ['bishop_cages_triangulated_embedding.msh']
outFiles = ['chessBishop_deformed.obj']

def eval_runtime(meshFile, cageFile, cageDeformedFile, embeddingFile, outFile, coordsFlag):
    command = ['cageDeformation3D', '-m', meshFile, '-c', cageFile, '--cd', cageDeformedFile, '-e', embeddingFile, '-o', outFile, coordsFlag, '-v', '0', '-t']
    pipe = subprocess.Popen(command, stdout=subprocess.PIPE)
    log = pipe.communicate()[0].decode("utf-8")
    if (pipe.returncode != 0):
        print('cageDeformation3D failed for ' + str(command))
        print(log)
    return float(log)

def eval_runtimes_meshes():
    f = open('runtimes.csv', 'w')
    writer = csv.writer(f)
    for i in range(len(models)):
        meshFile = models[i]
        cageFile = cages[i]
        cageDeformedFile = cagesDeformed[i]
        embeddingFile = embeddings[i]
        outFile = outFiles[i]
        runtimes_row = [meshFile]
        for coordsFlag in coordinateFlags:
            print('Evaluate ' + meshFile + ' (...) ' + coordsFlag)
            runtimes_row.append(eval_runtime(meshFile, cageFile, cageDeformedFile, embeddingFile, outFile, coordsFlag))
        writer.writerow(runtimes_row)

def main():
    os.chdir('../models')
    eval_runtimes_meshes()


if __name__ == "__main__":
    main()
