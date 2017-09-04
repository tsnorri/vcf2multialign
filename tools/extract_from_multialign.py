#
# Copyright (c) 2017 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import os
import sys

# Read from file into a buffer.
# Used the idea from https://stackoverflow.com/a/26209275/856976
def chunks(fp, bufsize):
	while True:
		chunk = fp.read(bufsize)
		if not chunk:
			break
		
		yield chunk


def chars(fp, bufsize = 4096):
	for chunk in chunks(fp, bufsize):
		for char in chunk:
			yield char
			

def handle_file(input, input_is_reference, offset, length, output, found_offset = None):
	k = 0
	if input_is_reference:
		i = 0
		j = 0
		outputting = False
		for char in chars(input):
			if char != '-':
				if i == offset:
					if found_offset is not None:
						found_offset(j)
					outputting = True
				i += 1
				
				if outputting:
					output.write(char)
					k += 1
					if k == length:
						break
			j += 1
	else:
		input.seek(offset, 0)
		for char in chars(input):
			output.write(char)
			k += 1
			if k == length:
				break


if __name__ == "__main__":
	parser = argparse.ArgumentParser("Extract subsequences from vcf2multialign output.")
	parser.add_argument('--input', type = argparse.FileType('rU'), required = True)
	parser.add_argument("--input-is-reference", action = 'store_true', default = False)
	parser.add_argument('--offset', type = int, required = True)
	parser.add_argument('--length', type = int, required = True)
	args = parser.parse_args()

	if args.offset < 0:
		parser.error("Offset has to be non-negative.")
	if args.length <= 0:
		parser.error("Length must be positive.")

	handle_file(
		args.input,
		args.input_is_reference,
		args.offset,
		args.length,
		sys.stdout,
		lambda n: print("Found the requested substring at file offset %d" % n, file = sys.stderr)
	)
