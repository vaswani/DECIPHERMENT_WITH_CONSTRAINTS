'''
converting from ARPA to carmel format'
'''
from optparse import OptionParser
import sys
import math
from collections import *
import re

order_re = re.compile(r'\\(\d)-grams:')

def main():
  # setting up options
  parser = OptionParser()
  parser.add_option("--ARPA", action="store", type="string", dest="arpa_file",default="lm.arpa",help="The LM in ARPA format")
  parser.add_option("--CARMEL", action="store", type="string", dest="carmel_file",default = 0.0,help="the carmel output file")
  parser.add_option("--order", action="store", type="int", dest="order",default = 3,help="the carmel output file")

  
  (options, args) = parser.parse_args()
  arpa_file = options.arpa_file
  carmel_file = options.carmel_file
  order = options.order
  one_grams_f = defaultdict(float)
  one_grams_bow = defaultdict(float)
  one_grams_p = defaultdict(float)
  two_grams_f = defaultdict(lambda :defaultdict(float))
  two_grams_bow = defaultdict(lambda :defaultdict(float))
  two_grams_p = defaultdict(lambda :defaultdict(float))
  three_grams_f = defaultdict(lambda :defaultdict(lambda :defaultdict(float)))
  three_grams_bow = defaultdict(lambda :defaultdict(lambda: defaultdict(float)))
  three_grams_p = defaultdict(lambda :defaultdict(lambda: defaultdict(float)))
  four_grams_p = defaultdict()
  ngrams_p = []
  ngrams_bow = []
  ngrams_f = []
  # lambda function for creating defaultdicts of a specified size
  g = lambda a:defaultdict(float) if a ==1 else defaultdict(lambda :g(a-1))
  for i in range(1,order+1) :
    print 'i is '
    print 'g(i)',g(i)
    ngrams_p.append(g(i))
    ngrams_f.append(g(i))
    ngrams_bow.append(g(i))
  '''
  print ngrams_p
  print ngrams_bow
  '''
  #raw_input()
  readArpaFile(arpa_file,ngrams_f,ngrams_bow,order)
  print 'read arpa file'
  #raw_input()
  #convertToProbabilities(one_grams_f,one_grams_bow,one_grams_p,two_grams_f,two_grams_bow,two_grams_p,\
  #    three_grams_f,three_grams_bow,three_grams_p)
  all_letters = []
  for k in range(65, 91):
    plain_letter = chr(k)
    all_letters.append(plain_letter)
  all_letters.append('_')
  print 'all letters is '
  print all_letters
  raw_input()

  writeFSA(carmel_file,ngrams_f,ngrams_bow,order,all_letters)
  ##a (aq *e* "q" 0.000426774195215!))

def writeFSA(carmel_file,ngrams_f,ngrams_bow,order,all_letters):
  output_file = open(carmel_file,'w')
  output_file.write("END\n")
  #for next_letter in two_grams_f['<s>']:
  #  output_file.write("(0 (#%s *e* \"%s\" %f!)))\n"%(next_letter,next_letter,math.pow(10,two_grams_f['<s>'][next_letter])))a
  two_grams_f = ngrams_f[1]
  for letter in two_grams_f['<s>']:
    if letter == '</s>':
      print 'ERROR! There was an empty path'
      sys.exit(0)
    else :
      if order == 3:
        output_file.write("(0 (#%s *e* \"%s\" %f!))\n"%(letter,letter,math.pow(10,two_grams_f['<s>'][letter])))
      elif order == 4:
        output_file.write("(0 (##%s *e* \"%s\" %f!))\n"%(letter,letter,math.pow(10,two_grams_f['<s>'][letter])))
  
  if order == 3 :
    three_grams_f = ngrams_f[2]
    for letter in three_grams_f:
      initial_state = []
      for next_letter in three_grams_f[letter]:
        initial_state = ''
        if letter == '<s>':
          initial_state = "#"
        else :
          initial_state = letter
        initial_state += next_letter
        for final_letter in three_grams_f[letter][next_letter] :
          final_state = ''
          if final_letter == '</s>':
            final_state= 'END'
            output_file.write("(%s (%s *e* *e*  %f!))\n"%(initial_state,final_state,math.pow(10,three_grams_f[letter][next_letter][final_letter])))
          else :
            final_state = "%s%s"%(next_letter,final_letter)
            output_file.write("(%s (%s *e* \"%s\" %f!))\n"%(initial_state,final_state,final_letter,math.pow(10,three_grams_f[letter][next_letter][final_letter])))
  elif order == 4:
    three_grams_f = ngrams_f[2]
    for letter in three_grams_f:
      initial_state = []
      for next_letter in three_grams_f[letter]:
        initial_state = ''
        if letter == '<s>':
          initial_state = "##"
        else :
          initial_state = "#%s"%letter
        initial_state += next_letter
        for final_letter in three_grams_f[letter][next_letter] :
          final_state = ''
          if final_letter == '</s>':
            final_state= 'END'
            output_file.write("(%s (%s *e* *e*  %f!))\n"%(initial_state,final_state,math.pow(10,three_grams_f[letter][next_letter][final_letter])))
          else :
            if letter == '<s>':
              final_state = "#%s%s"%(next_letter,final_letter)
            else :
              final_state = "%s%s%s"%(letter,next_letter,final_letter)
            output_file.write("(%s (%s *e* \"%s\" %f!))\n"%(initial_state,final_state,final_letter,math.pow(10,three_grams_f[letter][next_letter][final_letter])))
    four_grams_f = ngrams_f[3]
    for letter in four_grams_f:
      initial_state = []
      for next_letter_1 in four_grams_f[letter]:
        for next_letter_2 in four_grams_f[next_letter_1]:
          initial_state = ''
          if letter == '<s>':
            initial_state = "#"
          else :
            initial_state = "%s"%letter
          initial_state += next_letter_1
          initial_state += next_letter_2
          for final_letter in four_grams_f[letter][next_letter_1][next_letter_2] :
            final_state = ''
            if final_letter == '</s>':
              final_state= 'END'
              output_file.write("(%s (%s *e* *e*  %f!))\n"%(initial_state,final_state,math.pow(10,four_grams_f[letter][next_letter_1][next_letter_2][final_letter])))
            else :
              final_state = "%s%s%s"%(next_letter_1,next_letter_2,final_letter)
              output_file.write("(%s (%s *e* \"%s\" %f!))\n"%(initial_state,final_state,final_letter,math.pow(10,four_grams_f[letter][next_letter_1][next_letter_2][final_letter])))

  output_file.close()


  
def readArpaFile(arpa_file,ngrams_f,ngrams_bow,order):
  data_letters = {}
  f_dict = ngrams_f[0]
  bow_dict = ngrams_bow[0]
  start_reading = False
  for line in open(arpa_file):
    line = line.strip()
    if line == '':
      continue
    if r'-grams:' in line:
      print line
      current_order = 1
      order_match = order_re.search(line)
      if order_match == None:
        print 'error. Line ',line,' should have matched regular expression'
      else :
        current_order = int(order_match.group(1))
      start_reading = True
      f_dict = ngrams_f[current_order-1]
      bow_dict = ngrams_bow[current_order-1]
      print 'reading ',current_order,' grams'
      #raw_input()
      '''
      elif line == r'\2-grams:':
        f_dict = two_grams_f
        bow_dict = two_grams_bow
        print 'reading 2 grams'
        #raw_input()
      elif line == r'\3-grams:':
        f_dict = three_grams_f
        bow_dict = three_grams_bow
        print 'reading 3 grams'
        #raw_input()
      '''
    elif line == '\end\\' :
      continue
    elif start_reading == True :
      items = line.split('\t')
      #print items
      letters = items[1].split()
      populateDict(f_dict,letters,float(items[0]))
      if len(items) == 3:
        populateDict(bow_dict,letters,float(items[2]))
      for letter in letters:
        data_letters[letter] = 1
  print 'size of letter is ',len(data_letters)
  print data_letters
  raw_input()
  '''
  for i in range(order):
    print 'printing',i,'gram'
    raw_input()
    print ngrams_f[i]
    raw_input()
  '''
  '''
  print one_grams_f
  raw_input()
  print one_grams_bow
  raw_input()
  print two_grams_f
  raw_input()
  print two_grams_bow
  raw_input()
  print three_grams_f
  raw_input()
  print three_grams_bow
  raw_input()
  '''

  #for letter_a in letters:
  #  for letter_b in letters
  #now check for sum to one
  '''
  for letter in three_grams_f:
      for next_letter in three_grams_f[letter] :
        sum = 0.
        print 'number of last letters is ',len(three_grams_f[letter][next_letter])
        for last_letter in three_grams_f[letter][next_letter]:
          print' log prob 10 is ',three_grams_f[letter][next_letter][last_letter]
          sum += math.pow(10,three_grams_f[letter][next_letter][last_letter])
        print 'sum is ',sum
        raw_input()
  '''
def populateDict(arg_dict,letters,number) :
  if len(letters) == 1:
    arg_dict[letters[0]] = number
  elif len(letters) == 2 :
    arg_dict[letters[0]][letters[1]] = number
  elif len(letters) == 3 :
    arg_dict[letters[0]][letters[1]][letters[2]] = number
  elif len(letters) == 4 :
    arg_dict[letters[0]][letters[1]][letters[2]][letters[3]]= number


def convertToProbabilities(one_grams_f,one_grams_bow,one_grams_p,two_grams_f,two_grams_bow, two_grams_p,\
      three_grams_f,three_grams_bow,three_grams_p):
  for letter in two_grams_f:
    sum =0.
    for next_letter in two_grams_f :
      two_grams_p[letter][next_letter] = math.pow(10,two_grams_f[letter][next_letter]) - math.pow(10,one_grams_bow[next_letter]+one_grams_f[next_letter])
      sum += two_grams_p[letter][next_letter]
      print two_grams_p[letter][next_letter]
    print 'sum for ',letter,next_letter,' is ',two_grams_p[letter][next_letter]
    raw_input()

  
if __name__=="__main__":
  main()




