ORFcompare
===============

**ORFcompare** is a companion tool to ORFanage that compares CDS annotations between query and template transcripts.
While ORFanage is designed to annotate new ORFs, ORFcompare evaluates and quantifies the differences between 
existing CDS annotations in two GTF/GFF files.

This tool is particularly useful for:

1. Validating ORFanage results against known annotations
2. Comparing CDS annotations from different sources
3. Quantifying frame preservation between transcript isoforms
4. Identifying transcripts with matching or divergent coding regions

Basic Usage
----------------------------------------------------------

At minimum, **ORFcompare** requires a query file, a template file, and an output file:

::

	$ orfcompare --query query.gtf --template template.gtf --output comparison.tsv

To include start and stop codon information, provide a reference genome:

::

	$ orfcompare --reference genome.fa --query query.gtf --template template.gtf --output comparison.tsv


Output Format
-------------------------------------

**ORFcompare** outputs a tab-separated file containing metrics for each query/template pair comparison.
For each query transcript that overlaps with a template transcript, the following metrics are computed:

- **Length metrics**: CDS lengths for both query and template
- **Overlap metrics**: Matching, in-frame, out-of-frame, extra, and missing bases
- **Percent identity metrics**: LPI, MLPI, and ILPI scores
- **Codon information**: Start and stop codons for both transcripts (requires ``--reference``)

For detailed column descriptions, please refer to the :ref:`ORFcompare Stats Output <orfcompare-stats-file>` section in File Formats.


Description of Options
--------------------------------------

""""""""""""""""""""""""
``--query``
""""""""""""""""""""""""

Path to the GTF/GFF file containing query transcripts with CDS annotations to be compared.

""""""""""""""""""""""""
``--template``
""""""""""""""""""""""""

Path to the GTF/GFF file containing template (reference) transcripts with CDS annotations.

""""""""""""""""""""""""
``--output``
""""""""""""""""""""""""

Path to the output TSV file where comparison results will be written.

""""""""""""""""""""""""
``--reference``
""""""""""""""""""""""""

Path to the reference genome in FASTA format. When provided, enables extraction of actual 
start and stop codon amino acids for both query and template transcripts. This allows 
verification of proper translation initiation (M for methionine) and termination (``.`` for stop codon).

""""""""""""""""""""""""
``--threads``
""""""""""""""""""""""""

Number of threads to use for parallel processing. Similar to ORFanage, transcripts are grouped 
by coordinate overlap or gene ID and processed independently.

""""""""""""""""""""""""
``--use_id``
""""""""""""""""""""""""

When enabled, transcripts are grouped by gene ID rather than coordinate overlap. 
This is useful when gene IDs are consistent between query and template files.


Interpreting Results
--------------------------------------

The key metrics to evaluate CDS similarity are:

**ILPI (In-frame Length Percent Identity)**
    The most important metric for coding sequence comparison. High ILPI (>90%) indicates 
    the query and template share most of their coding sequence in the same reading frame.

**Start/Stop Codon Match**
    - ``M`` for start codon indicates canonical translation initiation
    - ``.`` for stop codon indicates proper termination
    - Other amino acid letters indicate incomplete or alternative codons

**len_extra / len_missing**
    These values quantify how the query CDS differs from the template:
    
    - High ``len_extra``: Query has additional coding sequence not in template
    - High ``len_missing``: Query is missing coding sequence present in template

Example Analysis
--------------------------------------

Compare ORFanage output against the original reference annotation:

::

	$ orfcompare --reference genome.fa \
	             --query orfanage_output.gtf \
	             --template original_annotation.gtf \
	             --output validation.tsv

Then filter for high-confidence matches:

.. code-block:: bash

	# Find transcripts with >95% in-frame identity
	awk -F'\t' '$12 > 95' validation.tsv > high_confidence.tsv
	
	# Find transcripts with matching start codons
	awk -F'\t' '$13 == "M" && $14 == "M"' validation.tsv > matching_starts.tsv

