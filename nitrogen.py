from load import load_nitrogenase_seq
from load import load_metagenome

def longestSubstrings(meta, nitrogen):
	ans = []
	for m in meta:
		length = 0
		longest = ''
		for x in range(len(nitrogen)):
			span = 1
			while(nitrogen[x:x+span] in m[1] and x+span < len(nitrogen)):
				if(span>length):
					length = span
					longest = nitrogen[x:x+span]
				span += 1
		if(length >= 20):
			ans.append([m[0], length, longest])
	return ans



if __name__ == '__main__':
	nitrogen = load_nitrogenase_seq().replace('\n', '')
	meta = load_metagenome()
	print(longestSubstrings(meta, nitrogen))