import sys

tag_file = open(sys.argv[1])

tags = []
for line in tag_file :
	line = line.strip()
	line_tags = line.split(' ')
	tags.append(line_tags)

tag_file.close()

bigrams = {}

for i in range(0,len(tags)) :
	for j in range(1,len(tags[i])) :
		if tags[i][j-1] in bigrams :
			if tags[i][j] in bigrams[tags[i][j-1]] :
				bigrams[tags[i][j-1]][tags[i][j]] += 1
			else :
				bigrams[tags[i][j-1]][tags[i][j]] = 1
		else :
			temp_dict = {}
			temp_dict[tags[i][j]] = 1
			bigrams[tags[i][j-1]] = temp_dict

for parent_tag in bigrams.keys() :
	total_count = 0
	for child_tag in bigrams[parent_tag].keys() :
		total_count += bigrams[parent_tag][child_tag]
	for child_tag in bigrams[parent_tag].keys() :
		bigrams[parent_tag][child_tag] = float(bigrams[parent_tag][child_tag])/float(total_count)


for parent_tag in bigrams.keys() :
	total_count = 0
	for child_tag in bigrams[parent_tag].keys() :
		print 'the language probability of  %s|%s  is %s'%(child_tag,parent_tag,str(bigrams[parent_tag][child_tag]))
