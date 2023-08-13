import os
import argparse
import pandas as pd
import genomenotebook as gn
from collections import defaultdict
from bokeh.models import Title


def parse_args():
    parser = argparse.ArgumentParser(usage='plot_operons.py -r REFERENCE -g GFF -a ANTIGENS -o OUTPUT',
                                     description='''Searches for the reference genes that are absent in the found operons.''')

    parser.add_argument('-g', '--gff', default=None, nargs=1,
                        help='Operons genes gff file')
    parser.add_argument('-a', '--antigens', default=None, nargs=1,
					    help='O-antiigens operons table')
    parser.add_argument('-o', '--output', default='operons_imgs', nargs=1,
                        help='Output directory')
    return parser.parse_args()
    
    
if __name__ == '__main__':
	gff_file = parse_args().gff[0]
	antigens_file = parse_args().antigens[0]
	output_dir = parse_args().output[0]
	
	if not os.path.isfile(gff_file):
		print('GFF file not found')
		sys.exit(1)
	if not os.path.isfile(antigens_file):
		print('O-antigens file not found')
		sys.exit(1)
       
	if os.path.isdir(output_dir):
		if os.listdir(output_dir):
			print(f'Warning. {output_dir} is not empty!')
	else:
		print(f'Creating {output_dir} directory')
		os.makedirs(output_dir)
   
	gff = gn.parse_gff(gff_file)
	operons = pd.read_csv(antigens_file, sep='\t')
	
	glyphs = gn.get_default_glyphs()
	glyphs = defaultdict()
	glyphs["CDS"] = gn.Glyph(glyph_type="arrow", 
		                    colors="blue", 
		                    height=0.8, 
		                    show_name=True)
	glyphs["transposon"] = gn.Glyph(glyph_type="arrow", 
		                    colors="red", 
		                    height=0.8, 
		                    show_name=True)
	
	for index, row in operons.iterrows():
		operon_genes_gff = gff[gff.operon == str(row.operon)]
		n_genes = len(operon_genes_gff)
		max_gene_name_len = operon_genes_gff.locus_tag.str.len().max()
		print(row.operon)
		g=gn.GenomeBrowser(gff_path=gff_file, #genome_path=genome_path,
		               glyphs=glyphs, 
		               search = False,
		               init_win = row.right + 3000 - row.left, 
		               bounds = (row.left - 1000, row.right + 1000),
		               feature_types= ['CDS', 'isertion_sequence'],
		               height=50 + max_gene_name_len*13 + n_genes * 26,
		               output_backend="svg",
		               feature_height=0.1,
		               attributes = ["locus_tag","product","start","end"], #will be displayed when hovering,
		               title=f"Operon {row.operon}"
		               )
		g.highlight(data=operons, hover_data=['operon', 'N_genes'])
		g.gene_track.title.align = 'center'
		g.gene_track.title.text_font_size = "25px"

		for gene_i, gene_row in operon_genes_gff.iterrows():
		    text = f'{gene_row.locus_tag}:  {gene_row["product"]}'
		    g.gene_track.add_layout(Title(text=text, text_font_size="9pt", text_font_style="normal"), 'below')
		    
		g.save(os.path.join(output_dir, f'{row.operon}.svg'))
		                    

	
	
