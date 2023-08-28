# Copyright (c) 2020 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import sys

inputs = [
	("test-1a.vcf",	"test-1.fa"),
	("test-1b.vcf",	"test-1.fa"),
	("test-2.vcf",	"test-2.fa"),
	("test-3.vcf",	"test-3.fa"),
	("test-4.vcf",	"test-4.fa")
]

for vcf, fasta in inputs:
	base = vcf.rstrip(".vcf")
	graph_name = f"{base}.graph"
	dot_name = f"{base}.dot"
	pdf_name = f"{base}.pdf"
	sys.stdout.write(f"../../../create-variant-graph/create_variant_graph -r {fasta} -a {vcf} -o {graph_name} -O -c 1\n")
	sys.stdout.write(f"../../../variant-graph-to-gv/variant_graph_to_gv -r {fasta} -g {graph_name} -o {dot_name}\n")
	sys.stdout.write(f"dot -Tpdf -o{pdf_name} {dot_name}\n")
