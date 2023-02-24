
#***************************************************************************************************************************************#
#                                                                                                                                       #
#       Script name:    task_2.py                                                                                                       #
#       Description:    The NCBI database records also sequences for SARS-COV: https:// www.ncbi.nlm.nih.gov/sars-cov-2/ 
#                       This task is about using the NCBI databases and BLAST to compare the amino acid sequence of SARS-CoV-2, 
#                       SARS-CoV and MERS-CoV.                                                          #                     
#       Module:         Algorithms and Data Structures in Bioinformatics - ADS 612                                                      #
#       Submitted By:   Vita Nyasulu(202240190002), Samir Adrian Harawa(202240190001) & Limbani Thengo(202240190021)                                                                     #
#       Date:           24 February 2023                                                                                                #
#                                                                                                                                       #
#***************************************************************************************************************************************#


# Import necessary libraries
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
import ssl

# Setting up the context for secure connection
ssl._create_default_https_context = ssl._create_unverified_context

# Read the SARS-CoV-2 genome file
record = SeqIO.read("wuhan.fasta", "fasta")

# Extract the nucleotide and amino acid sequences for SARS-CoV-2, SARS-CoV, and MERS-CoV from the NCBI database
# Set the email address to use when making requests to NCBI

Entrez.email = "vitanyasulu@gmail.com" # replace with your email address
# Retrieve the protein sequence for SARS-CoV-2 from the NCBI protein database using its ID
handle_sars2 = Entrez.efetch(db="protein", id="YP_009724390.1", rettype="fasta", retmode="text")
record_sars2 = SeqIO.read(handle_sars2, "fasta")
seq_sars2 = str(record_sars2.seq)
# Retrieve the protein sequence for SARS-CoV from the NCBI protein database using its ID
handle_sars = Entrez.efetch(db="protein", id="AAP13441.1", rettype="fasta", retmode="text")
record_sars = SeqIO.read(handle_sars, "fasta")
seq_sars = str(record_sars.seq)
# Retrieve the protein sequence for MERS-CoV from the NCBI protein database using its ID
handle_mers = Entrez.efetch(db="protein", id="YP_009047204.1", rettype="fasta", retmode="text")
record_mers = SeqIO.read(handle_mers, "fasta")
seq_mers = str(record_mers.seq)

# Run a BLAST search to find the most similar protein sequences in the NCBI database
result_handle = NCBIWWW.qblast("blastp", "nr", seq_sars2)
blast_record = NCBIXML.read(result_handle)

# Print the top 10 hits
for alignment in blast_record.alignments[:10]:
    # Print information about each hit alignment
    print("****Alignment****")
    print("sequence title:", alignment.title[:80])
    print("sequence length:", alignment.length)
    for hsp in alignment.hsps:
        # Print information about each hit alignment's high-scoring segment pair (HSP)
        print("e value:", hsp.expect)
        print(hsp.query[:80] + "...")
        print(hsp.match[:80] + "...")
        print(hsp.sbjct[:80] + "...")


# Use the Biopython GenomeDiagram to visualize the ORFs in the SARS-CoV-2 genome
# Create a new diagram object for the SARS-CoV-2 genome
gd_diagram = GenomeDiagram.Diagram("SARS-CoV-2")
#Create a new track for the annotated features and add it to the diagram
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
# Create a new feature set for the annotated features track
gd_feature_set = gd_track_for_features.new_set()
# Iterate over each feature in the SARS-CoV-2 genome
for feature in record.features:
    # Check if the feature is a CDS (coding sequence)
    if feature.type == "CDS":
        # Set the color for the feature
        color = colors.lightblue
        # Get the name of the protein from the feature's qualifiers or set it as "unknown
        name = feature.qualifiers.get("product", ["unknown protein"])[0]
        gd_feature = SeqFeature(FeatureLocation(feature.location.start, feature.location.end), strand=feature.location.strand)
        gd_feature_set.add_feature(gd_feature, name=name, color=color)
# Save the GenomeDiagram as a PNG file
gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                start=0, end=len(record), circle_core=0.5)

gd_diagram.write("SARS-CoV-2_ORFs.png", "PNG")
