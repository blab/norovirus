import os
import sys
import fileinput
# Read in the file
with open('metadata_parsed.tsv', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('Viet_Nam', 'Vietnam')

# Write the file out again
with open('metadata_parsed.tsv', 'w') as file:
  file.write(filedata)
