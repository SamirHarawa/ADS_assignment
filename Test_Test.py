from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import BiopythonTranslator
import ssl

# Setting up the context for secure connection
ssl._create_default_https_context = ssl._create_unverified_context

# Read the SARS-CoV-2 genome file
record = SeqIO.read("wuhan.fasta", "fasta")

# Calculate the GC% content
gc_count = record.seq.count("G") + record.seq.count("C")
gc_percent = (float(gc_count) / len(record.seq)) * 100

# Generate a graph showing the GC% content
gc_graph = GenomeDiagram.Diagram("GC Content")
gc_track = gc_graph.new_track(1, name="GC Content")
gc_feature_set = gc_track.new_set()

gc_feature = SeqFeature(FeatureLocation(0, len(record)), strand=None)
gc_feature_set.add_feature(gc_feature, color=colors.green, name="GC Content")

gc_graph.draw(format="linear", pagesize=(20 * cm, 5 * cm))
gc_graph.write("gc_content.png", "PNG")

# Extract the nucleotide and amino acid sequences for SARS-CoV-2, SARS-CoV, and MERS-CoV from the NCBI database
Entrez.email = "vitanyasulu@gmail.com"  # replace with your email address
handle_sars2 = Entrez.efetch(db="protein", id="YP_009724390.1", rettype="fasta", retmode="text")
record_sars2 = SeqIO.read(handle_sars2, "fasta")
seq_sars2 = str(record_sars2.seq)

handle_sars = Entrez.efetch(db="protein", id="AAP13441.1", rettype="fasta", retmode="text")
record_sars = SeqIO.read(handle_sars, "fasta")
seq_sars = str(record_sars.seq)

handle_mers = Entrez.efetch(db="protein", id="YP_009047204.1", rettype="fasta", retmode="text")
record_mers = SeqIO.read(handle_mers, "fasta")
seq_mers = str(record_mers.seq)

# Run a BLAST search to find the most similar protein sequences in the NCBI database
result_handle = NCBIWWW.qblast("blastp", "nr", seq_sars2)
blast_record = NCBIXML.read(result_handle)

# Print the top 10 hits
for alignment in blast_record.alignments[:1]:
    print("****Alignment****")
    print("sequence title:", alignment.title[:80])
    print("sequence length:", alignment.length)
    for hsp in alignment.hsps:
        print("e value:", hsp.expect)
        print(hsp.query[:80] + "...")
        print(hsp.match[:80] + "...")
        print(hsp.sbjct[:80] + "...")

# Use the dna_features_viewer library to visualize the the ORFs in the SARS-CoV-2 genome
from Bio.Graphics import GenomeDiagram
from dna_features_viewer import GraphicFeature, CircularGraphicRecord

gd_diagram = BiopythonTranslator().translate_record(record)

# Create a diagram object
diagram = GenomeDiagram.Diagram("SARS-CoV-2_ORFs")

# Add the track and the features to the diagram
features=[
    GraphicFeature(start=0, end=20, strand=+1, color="#ffd700",
                   label="Small feature"),
    GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a very long name"),
    GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                   label="Gene 2"),
    GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
                   label="Gene 3")
]
circ_record = CircularGraphicRecord(sequence_length=len(record), features=features)
circ_record.plot(figure_width=5)

gc_track = gc_graph.new_track(1, name="GC Content")
diagram.add_track(gc_track, 1)
feature_set = gc_track.new_set()

for feature in gd_diagram.features:
    if feature.type == "gene":
        if feature.strand > 0:
            color = colors.blue
        else:
            color = colors.red
        feature_set.add_feature(feature, color=color, label=True)

# Draw the diagram with a circular layout
diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm), start=0, end=len(record), circle_core=0.5,track_size=0.4)
diagram.write("SARS-CoV-2_ORFs.png", "PNG")

#gd_diagram.plot("SARS-CoV-2_ORFs.png")

# Find the coding regions of CDS in the Covid-19 genome
cds_features = [f for f in record.features if f.type == 'CDS']

# ORFs identification, find the Open Reading Frames in the COVID-19 genome, set the minimum protein length to 200 amino acids
orf_features = []
for f in cds_features:
    if len(f.qualifiers.get("translation"))[0] >= 200:
        orf_features.append(f)

#Add the spike protein sequence to the graphic
spike_seq = seq_sars2[21563:25384]
spike_feature = SeqFeature(FeatureLocation(21563, 25384), strand=None)
orf_features.append(spike_feature)

#Add features to the graphic
gd_diagram_features = []
for feature in orf_features:
    if feature.location.strand == 1:
        color = colors.blue
    else:
        color = colors.red
        gd_feature = feature.extract(record)
        gd_feature.color = color
        gd_diagram_features.append(gd_feature)

#Add a feature for the spike protein
spike_gd_feature = SeqFeature(FeatureLocation(21563, 25384), strand=None)
spike_gd_feature.color = colors.red
gd_diagram_features.append(spike_gd_feature)

#Add features for the S, E, M, and N proteins
s_feature = SeqFeature(FeatureLocation(265, 13468), strand=None)
s_feature.color = colors.blue
gd_diagram_features.append(s_feature)

e_feature = SeqFeature(FeatureLocation(26244, 26472), strand=None)
e_feature.color = colors.purple
gd_diagram_features.append(e_feature)

m_feature = SeqFeature(FeatureLocation(26523, 27191), strand=None)
m_feature.color = colors.green
gd_diagram_features.append(m_feature)

n_feature = SeqFeature(FeatureLocation(28274, 29533), strand=None)
n_feature.color = colors.yellow
gd_diagram_features.append(n_feature)

#Render the diagram
gd_diagram = BiopythonTranslator().translate_record(record)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

for feature in gd_diagram.features:
    if feature.type == "gene":
        if feature.strand > 0:
            color = colors.blue
        else:
            color = colors.red
        gd_feature_set.add_feature(feature, color=color, label=True)

#Frequency of the Amino acid in the spike protein
aa_count = {}
for aa in spike_seq:
    if aa in aa_count:
        aa_count[aa] += 1
    else:
        aa_count[aa] = 1

#Calculate the frequency of each amino acid
total_aa = sum(aa_count.values())
aa_freq = {k: (v / total_aa) for k, v in aa_count.items()}

#Print the frequency of each
for aa, freq in aa_freq.items():
    print(f"{aa}: {freq:.2%}")

#Plot the frequency of each amino acid in the spike protein
import matplotlib.pyplot as plt
plt.bar(aa_freq.keys(), aa_freq.values())
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.title("Spike Protein Amino Acid Frequency")
plt.show()

#Calculate the total number of amino acids in the spike protein
total_aa_spike = sum(aa_count.values())
print("Total number of amino acids in the spike protein:", total_aa_spike)