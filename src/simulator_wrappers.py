import sys
import os

class Simulator:
    "art_illumina -sam -i reference.fa -p -l 150 -ss HS25 -f 20 -m 200 -s 10 -o paired_dat"
    
    @staticmethod
    def get_art_illumina(genome, output_prefix, paired_end, read_length, vertical_coverage, mean_fragment_size=200, stddev_insert=10):
        result = ["art_illumina"]
        if paired_end:
            result.append("-p")
        result.append("-i")
        result.append(genome)
        result.append("-o")
        result.append(output_prefix)
        result.append("-l")
        result.append(read_length)
        result.append("-ss")
        result.append("HS25")
        result.append("-f")
        result.append(vertical_coverage)
        result.append("-m")
        result.append(mean_fragment_size)
        result.append("-s")
        result.append(stddev_insert)
        
        read1 = "{}_1.fq".format(output_prefix)
        read2 = "{}_2.fq".format(output_prefix)

        return result, read1, read2
    