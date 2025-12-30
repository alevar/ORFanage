.. _file-formats:

File Formats
======================

**ORFanage** is designed to use simple and common file formats.

.. _gtf-file:

Re-Annotated GTF Output
-------------------------

.. csv-table:: Re-annotated GTF Output (TAB-separated)
   :file: ../content/csvs/gtf.t.csv
   :header-rows: 0
   :stub-columns: 1


.. csv-table:: Attributes (column #9) (TAB-separated)
   :file: ../content/csvs/attributes.csv
   :header-rows: 0
   :stub-columns: 1


.. _stats-file:

ORFanage Stats Output (TSV)
---------------------------

The stats output file (``--stats``) contains detailed metrics for each query/template comparison.

.. csv-table:: ORFanage Stats Output (TSV)
   :file: ../content/csvs/stats.t.csv
   :header-rows: 1
   :stub-columns: 1

**Key Metrics:**

- **length_pi (LPI)**: Measures how much of the combined query+template CDS region is covered by the query. Formula: ``100 * query_len / union_len``.
- **match_length_pi (MLPI)**: Measures actual overlap between query and template, excluding query extensions. Formula: ``100 * len_match / union_len``.
- **inframe_length_pi (ILPI)**: The most biologically meaningful metric—measures what fraction maintains the correct reading frame. Formula: ``100 * len_inframe / union_len``.
- **pi**: Requires genome sequence (``--reference``); measures actual nucleotide sequence identity in the aligned regions at the codon level.


.. _orfcompare-stats-file:

ORFcompare Stats Output (TSV)
-----------------------------

The orfcompare output file compares CDS annotations between query and template transcripts.

.. csv-table:: ORFcompare Stats Output (TSV)
   :file: ../content/csvs/orfcompare_stats.csv
   :header-rows: 1
   :stub-columns: 1

**Codon Columns (require --reference):**

- **query_start_codon / template_start_codon**: Amino acid at the CDS start. ``M`` (methionine) indicates a canonical start codon.
- **query_stop_codon / template_stop_codon**: Amino acid at the CDS end. ``.`` indicates a proper stop codon; any other letter indicates the ORF is incomplete.
