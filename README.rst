==========================================
Carrying out GO term enrichment with topGO
==========================================

I'm not a big fan of R, but this is the one thing that I have not been able to replicate in Python--how to carry out a GO term enrichment that is more complicated than Fisher's exact testing + multiple-testing correction. Sigh.

As such, while my R code is *functional*, do not use it to learn how to code better in R. It's far... uglier than my usual standards.

This pipeline requires:

1. Read (and cite, if you find it useful) the paper to get an idea of how things work: https://academic.oup.com/bioinformatics/article/22/13/1600/193669.

2. R (the newer, the better. From experience, topGO has been continually updated to work with the newest R versions).

   Optional: install RStudio. Makes programming in R... bearable.

3. Install topGO via Bioconductor: http://bioconductor.org/packages/release/bioc/html/topGO.html.

4. Troubleshooting if things fail!

   Vignettes: http://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

   Manual: http://bioconductor.org/packages/release/bioc/manuals/topGO/man/topGO.pdf

Pipeline summary
----------------
1. Test whether topGO has been correctly installed.
2. Set up a folder to contain the necessary files.
3. Understanding when to pay attention to the gene universe, and when not to.
4. Run the scripts!

Test topGO installation
-----------------------
Go to R, then run ``library(topGO)``.

No errors? Yay!

Errors popped up? Google the error message, check whether your R is outdated, and reinstall the package.

Setting up relevant files and folders
-------------------------------------
In essence, this step is necessary because my R programming is terrible. You'll notice that the topGO R script contains a lot of hardcoded filenames and folders. That's mainly because I got lazy and didn't bother figuring out how to configure an R script to take in command-line arguments, like any sensible script should. It's not easy to pass arguments to the script via RStudio, hence my laziness.

OK, enough griping. I'll use an example (*Aiptasia*) to illustrate the process. If you mimic the way I set up this git repository, you'd be fine. Hopefully.

Create a folder (I like to call mine "go_annot").

Within "go_annot", there should be two folders and a bunch of loose files.

Loose files, in alphabetical order:

- ``aip_go_annots.all.tsv``: a plain-text, tab-separated file that contains gene names and all GO terms assigned to the respective genes.

  Sidenote: the file I provide here has a lot of GO terms, which might look strange to some of you. This is because in my annotation pipeline, if a gene is, for instance, assigned a level-10 GO term, my script recursively checks for the parent of that term, all the way till it hits one of the three level-1 GO terms ("biological process", "cellular component", "molecular function").
  
  Does this affect topGO? Nope! From testing, *almost exact* results are produced from using a fully-expanded GO term file, or if you used the most specific GO terms that applied to the gene. I believe it's because topGO, like what I did, recursively checks for parents of all terms, using the most recent GO term hierarchy. As to why I emphasise "almost exact", I think differences are introduced when the hierarchy of GO terms change over time--my expansions are frozen in time, but topGO's will continually get updated over time. Differences are really minor though, so don't lose sleep over it.
  
  If you don't really understand this bit, don't worry. Bottom line is, make sure your file has gene names and GO terms, and you're good to go.

- ``aip_topgo_usage.R``: the R script that you should modify (e.g. changing filenames) if you don't care about gene universes.

- ``aip_topgo_usage.consider_universe.R``: the R script you should use if you care about gene universes. I'm devoting a section later to explain what gene universes are; for now, if you're a first timer, ignore this file, use the previous one.

The two folders are:

1. ``genes_of_interest``

   This folder would contain subfolders containing your genes of interest, in plaintext. I have provided ``tx_vs_prot`` as an example, which contains four plaintext files within it.

2. ``topGO_output``

   This would be where the topGO output files are produced. For now, this should only contain ``summarise_topGO_output.sh``. Ignore the subfolder ``tx_vs_prot`` for now.

Understanding gene universes
----------------------------
Let's quickly go through the idea behind GO term enrichments. Imagine your organism as a bag of coloured M&Ms, and the genes as the M&Ms themselves. The entire bag of M&Ms is the **gene universe**.

Say you grab a bunch of M&Ms, 100 of them in total. You notice 51 of them red, 20 yellow, 20 green, 9 brown. Oh my god--51% reds! Before you get too excited, take a step back and think: is this surprising? It depends!

If your **gene universe** has, say, 50% reds, then no, it's not surprising to get 51 reds amongst your 100.

If your **gene universe** has, say, equal number of all colours (i.e. 25% reds), observing 51 in 100 is very unexpected.

Anecdotally, most of the papers I've read that perform GO term enrichments sets the **gene universe** to all genes in the organism. Is this a reasonable assumption? Again--it depends!

If they are looking at RNA-seq data, then yes, I agree that this is a reasonable assumption. The vast majority of genes are transcribed in any tested condition, so it is perfectly fine to approximate the **gene universe** as all genes in the organism.

If they are looking at, say, whole proteome sequencing, then nope, this assumption breaks down. It is much harder to extract all proteins in an organism. Most studies report the detection of peptides from ~10% of all genes. This means that the set of *detected genes* would have a different distribution of GO terms (= M&M colours) than *all genes*. Going deeper, this means that if we are looking at a subset of genes in the detected genes (e.g. differentially expressed proteins), the **gene universe** should be the set of detected genes because you can't observed differential protein expression of proteins that are not detected.

This was why I wrote ``aip_topgo_usage.consider_universe.R``.

Running the analysis
--------------------
1. Open ``aip_topgo_usage.R`` or ``aip_topgo_usage.consider_universe.R`` in plain ol' R, or RStudio (recommended).

2. Modify the folder names to suit your usage. These are in lines 2 and 6 of the scripts.

3. Run it!

4. It should take 10-20 minutes to run. Once it's done, you'll see a bunch of text files created in the folder ``topGO_output``. These files will start with ``bp_``, ``cc_`` and ``mf_``.

5. Create a folder (in my example, I called it ``tx_vs_prot``), then move all the text files into the folder. Leave the shell script outside.

6. To produce the ``summary_`` file, you have to use the command line to invoke the shell script. Navigate to the folder containing all the text files, then run

   ``../summarise_topGO_output.sh``

   This script implements the default cutoffs (p < 0.05, and the GO term has to exist at least 5 times in the **gene universe**). The summary files can then be used--terms in this file are ENRICHED (no, not depleted) in your gene of interest list).
   
   If this doesn't work for you, this can be done manually. Open the text files in Excel, then manually select for lines with the 4th column >= 5, and 7th column < 0.05.

Naming files for ``*.consider_universe.R``
------------------------------------------
Take a close look at how I named my gene of interest files in ``genes_of_interest/tx_vs_prot``. There are files called ``*_universe.txt``, and some called ``*_diff.txt``.

See how similar the filenames are? The ``*.consider_universe.R`` script knows that the **gene universe** for the ``cc7_prot_diff.txt`` is ``cc7_prot_universe.txt``, while ``cc7_RNA_diff.txt`` goes with ``cc7_RNA_universe.txt``, because it does a search-and-replace for "diff" to "universe". Yes, hacky, I know. In the script, you can find it at line 18--feel free to customise it to your liking. The in-built search-and-replaces are "up", "down" and "diff" to "universe".

Good luck!
