#
# Copyright (c) 2017 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import codecs
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


def handle_ref_input(src, offset, length):
	"""Find the offset of the given gene in the reference input."""
	src.seek(0, 0)
	i = 0
	j = 0
	outputting = False
	for char in chars(src):
		if char != '-':
			if i == offset:
				return j
			i += 1
		j += 1
	return None


def write_sequence(src, dst, fname, file_offset, length):
	k = 0
	src.seek(file_offset, 0)
	for char in chars(src):
		if '-' == char:
			continue

		dst.write(char)
		k += 1
		if k == length:
			break
			

def handle_files(ref_input, co_ordinate_input, seq_file_names):
	write_flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
	for line in co_ordinate_input:
		fields = line.strip().split("\t")
		fields = fields[0:6]
		(chrom, chrom_start, chrom_end, name, score, strand) = fields
		chrom_start = int(chrom_start)
		chrom_end = int(chrom_end)
		length = chrom_end - chrom_start
		print("Handling sequence '%s'…" % name, file = sys.stderr)

		# Find the offset from the reference file.
		file_offset = handle_ref_input(ref_input, chrom_start, length)
		print("Found the requested substring at file offset %d" % file_offset, file = sys.stderr)

		# Handle the source files.
		fd = os.open("%s.fa" % name, write_flags)
		with os.fdopen(fd, 'w') as dst:
			for fname in seq_file_names:
				with open(fname, 'r') as src:
					print("\tHandling source file '%s'…" % fname, file = sys.stderr)
					dst.write(">%s\n" % fname)
					write_sequence(src, dst, fname, file_offset, length)
					dst.write("\n")


if __name__ == "__main__":
	parser = argparse.ArgumentParser("Extract subsequences from vcf2multialign output.")
	parser.add_argument('--aligned-reference', type = argparse.FileType('rU'), required = True)
	parser.add_argument('--extracted-co-ordinates', type = argparse.FileType('rU'), required = True)
	parser.add_argument('source-files', nargs = '*')
	args = parser.parse_args()

	# Output UTF-8, https://stackoverflow.com/a/4374457/856976
	#sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())
	#sys.stderr = codecs.getwriter("utf-8")(sys.stderr.detach())

	handle_files(
		args.aligned_reference,
		args.extracted_co_ordinates,
		vars(args)['source-files'],
	)
