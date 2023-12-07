# Copyright (c) 2020-2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import sys

if __name__ == "__main__":
	inputs = [
		("test-1a.vcf",	"test-1.fa"),
		("test-1b.vcf",	"test-1.fa"),
		("test-2.vcf",	"test-2.fa"),
		("test-3.vcf",	"test-3.fa"),
		("test-4.vcf",	"test-4.fa")
	]
	
	for vcf, fasta in inputs:
		base = vcf.rstrip(".vcf")
		dot_name = f"{base}.dot"
		pdf_name = f"{base}.pdf"
		sys.stdout.write(f"../../../vcf2multialign/vcf2multialign -r {fasta} -a {vcf} -d {dot_name} -c 1\n")
		sys.stdout.write(f"dot -Tpdf -o{pdf_name} {dot_name}\n")
