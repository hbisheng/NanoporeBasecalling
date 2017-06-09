from subprocess import Popen, PIPE

#args = ["python /home/ubuntu/basecall/deepnano/basecall.py --output /home/ubuntu/dataplace/partial_data/deepnano7.3.fa --directory /home/ubuntu/dataplace/partial_data/R7.3/ > /home/ubuntu/dataplace/partial_data/deepnano7.3.log"]

#args = ["python", "basecall.py", "--directory", "/home/ubuntu/dataplace/samples7.3/"]
args = ["ls"]
p = Popen(['time', '-f', '%U user %S system %E elapsed %M Kbytes'] + args, stderr=PIPE)

print p.communicate()
