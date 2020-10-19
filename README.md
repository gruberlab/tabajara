# TABAJARA - Tool for Alignment Block Analysis Joining Appropriate Rational Approaches

TABAJARA is a tool for rational design of profile HMMs. Starting from a multiple sequence alignment (MSA), TABAJARA is able to find blocks that are either (1) conserved across all sequences or (2) discriminative for two specific groups of sequences. For the identification of regions conserved across all protein sequences of an MSA, we implemented a previously described algorithm (Capra & Singh, 2007), based on Jensen–Shannon divergence method (Lin, 1991). This is the method of choice to determine the level of character conservation across all sequences of an MSA. In the case of nucleotide sequences, TABAJARA uses Shannon entropy (Shannon, 1948; Shannon & Weaver, 1949) to calculate position-specific scores. To find group-discriminative blocks, the program can use either Mutual Information (Adami, 2004; Cover & Thomas, 2006) or Sequence Harmony (Feenstra et al., 2007; Pirovano et al., 2006) for both, DNA or protein sequences.
