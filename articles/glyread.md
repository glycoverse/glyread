# Get Started with glyread

## The Great Glycoproteomics Data Chaos 📊

Picture this: you’ve just finished running your glycoproteomics
experiment through your favorite identification software—maybe pGlyco3,
MSFragger-Glyco, or Byonic. You’re excited to dive into the results, but
then reality hits: you’re staring at a massive spreadsheet with 120+
columns, cryptic column names, and data scattered everywhere like
confetti after a celebration gone wrong!

Sound familiar? Welcome to the wonderful world of glycoproteomics data
formats!

Each software tool speaks its own unique “dialect”—what one calls
“Proteins,” another might call “Protein.Accessions,” and yet another
prefers “UniProt_IDs.” It’s like trying to have a conversation at the
Tower of Babel, but for data scientists.

**Enter glyread—your universal translator! 🌟**

Think of `glyread` as your personal data butler who speaks fluent
“software chaos” and translates everything into beautiful, organized,
analysis-ready data. No more wrestling with column names, no more manual
reformatting, no more pulling your hair out over inconsistent formats.
Just clean, tidy data that’s ready for the fun part: **actual
analysis!**

**🎯 Important Note:** All functions in `glyread` return a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object—the lingua franca of the `glycoverse` ecosystem. If you haven’t
met this elegant data structure yet, we highly recommend taking a quick
detour to [its
introduction](https://glycoverse.github.io/glyexp/articles/glyexp.html)
first. Trust us, it’s worth the journey! 🚀

``` r
library(glyread)
library(glyexp)
library(readr)
```

## A Real-World Example: Taming the pGlyco3 + pGlycoQuant Beast 🐉

Let’s dive into a classic glycoproteomics workflow that’ll make you
appreciate `glyread`’s magic!

**The Setup:** You’re using pGlyco3 (the glycopeptide identification
wizard 🧙‍♂️) paired with pGlycoQuant (the quantification maestro 🎼). This
dynamic duo is incredibly powerful, but their output? Well, let’s just
say it’s… comprehensive.

Here’s what a typical pGlycoQuant output file looks like:

``` r
read_tsv("glycopeptides.list")
#> New names:
#> Rows: 500 Columns: 120
#> ── Column specification
#> ──────────────────────────────────────────────────────── Delimiter: "\t" chr
#> (13): GlySpec, PepSpec, RawName, Peptide, Mod, Glycan(H,N,A,F), GlycanC... dbl
#> (105): Scan, RT, PrecursorMH, PrecursorMZ, Charge, Rank, PeptideMH, GlyI... lgl
#> (2): Empty_Separator, ...120
#> ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
#> Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> • `` -> `...120`
#> # A tibble: 500 × 120
#>    GlySpec      PepSpec RawName  Scan    RT PrecursorMH PrecursorMZ Charge  Rank
#>    <chr>        <chr>   <chr>   <dbl> <dbl>       <dbl>       <dbl>  <dbl> <dbl>
#>  1 20241224-LX… 202412… 202412…   266  58.1       2880.        961.      3     1
#>  2 20241224-LX… 202412… 202412…  2870 621.        3537.        885.      4     1
#>  3 20241224-LX… 202412… 202412…  2878 622.        3246.        812.      4     1
#>  4 20241224-LX… 202412… 202412…  2884 623.        3537.        708.      5     1
#>  5 20241224-LX… 202412… 202412…  3015 642.        2932.        978.      3     1
#>  6 20241224-LX… 202412… 202412…  3079 652.        3829.        958.      4     1
#>  7 20241224-LX… 202412… 202412…  3118 657.        3537.        885.      4     1
#>  8 20241224-LX… 202412… 202412…  3121 657.        3829.        767.      5     1
#>  9 20241224-LX… 202412… 202412…  3129 658.        3537.        708.      5     1
#> 10 20241224-LX… 202412… 202412…  3137 659.        3246.        812.      4     1
#> # ℹ 490 more rows
#> # ℹ 111 more variables: Peptide <chr>, Mod <chr>, PeptideMH <dbl>,
#> #   `Glycan(H,N,A,F)` <chr>, GlycanComposition <chr>, PlausibleStruct <chr>,
#> #   GlyID <dbl>, GlyFrag <chr>, GlyMass <dbl>, GlySite <dbl>, TotalScore <dbl>,
#> #   PepScore <dbl>, GlyScore <dbl>, CoreMatched <dbl>, MassDeviation <dbl>,
#> #   PPM <dbl>, GlyIonRatio <dbl>, byIonRatio <dbl>, czIonRatio <dbl>,
#> #   GlyDecoy <dbl>, PepDecoy <dbl>, Ion_163.06 <dbl>, Ion_366.14 <dbl>, …
```

That’s 120 columns of pure, unadulterated data chaos! While this
comprehensive output is fantastic for software interoperability, it’s
like trying to find a needle in a haystack when you just want to get to
your analysis.

**What you actually need** (the needle in our haystack):

1.  🧬 **Glycoform descriptions**: proteins, sites, glycan compositions,
    glycan structures
2.  📊 **Quantification results**: the actual numbers from pGlycoQuant

**The old way:** Manually wrestling with this 120-column monster,
writing custom parsing scripts, debugging column name mismatches, and
generally questioning your life choices.

**The glyread way:** One function call. Seriously. ✨

Meet
[`read_pglyco3_pglycoquant()`](https://glycoverse.github.io/glyread/reference/read_pglyco3_pglycoquant.md)—your
new best friend:

``` r
exp <- read_pglyco3_pglycoquant("glycopeptides.list", sample_info = "sample_info.csv")
#> ℹ Reading data
#> ℹ Finding leader proteins
#> ✔ Finding leader proteins [84ms]
#> 
#> ℹ Reading dataColumn group converted to <factor>.ℹ Parsing glycan compositions and structures
#> Column group converted to <factor>.✔ Parsing glycan compositions and structures [412ms]
#> 
#> ℹ Reading data✔ Reading data [824ms]
exp
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 298 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: peptide <chr>, peptide_site <int>, protein <chr>, protein_site <int>, gene <chr>, glycan_composition <comp>
```

**Ta-da!** Look at that beautiful transformation! From 120-column chaos
to organized elegance in one line of code.

Let’s peek at what treasures we’ve extracted:

**🏷️ Variable Information - Meet Your Glycopeptides:**

``` r
get_var_info(exp)
#> # A tibble: 298 × 7
#>    variable   peptide peptide_site protein protein_site gene  glycan_composition
#>    <chr>      <chr>          <int> <chr>          <int> <chr> <comp>            
#>  1 P08185-N1… NKTQGK             1 P08185           176 SERP… Hex(5)HexNAc(4)Ne…
#>  2 P04196-N3… HSHNNJ…            5 P04196           344 HRG   Hex(5)HexNAc(4)Ne…
#>  3 P04196-N3… HSHNNJ…            5 P04196           344 HRG   Hex(5)HexNAc(4)   
#>  4 P04196-N3… HSHNNJ…            5 P04196           344 HRG   Hex(5)HexNAc(4)Ne…
#>  5 P10909-N2… HNSTGC…            2 P10909           291 CLU   Hex(6)HexNAc(5)   
#>  6 P04196-N3… HSHNNJ…            5 P04196           344 HRG   Hex(5)HexNAc(4)Ne…
#>  7 P04196-J3… HSHNNJ…            6 P04196           345 HRG   Hex(5)HexNAc(4)   
#>  8 P04196-N3… HSHNNJ…            5 P04196           344 HRG   Hex(5)HexNAc(4)dH…
#>  9 P04196-N3… HSHNNJ…            5 P04196           344 HRG   Hex(4)HexNAc(3)   
#> 10 P04196-N3… HSHNNJ…            5 P04196           344 HRG   Hex(4)HexNAc(4)Ne…
#> # ℹ 288 more rows
```

**📋 Sample Information - Know Your Experiments:**

``` r
get_sample_info(exp)
#> # A tibble: 12 × 2
#>    sample                  group
#>    <chr>                   <fct>
#>  1 20241224-LXJ-Nglyco-H_1 H    
#>  2 20241224-LXJ-Nglyco-H_2 H    
#>  3 20241224-LXJ-Nglyco-H_3 H    
#>  4 20241224-LXJ-Nglyco-M_1 M    
#>  5 20241224-LXJ-Nglyco-M_2 M    
#>  6 20241224-LXJ-Nglyco-M_3 M    
#>  7 20241224-LXJ-Nglyco-Y_1 Y    
#>  8 20241224-LXJ-Nglyco-Y_2 Y    
#>  9 20241224-LXJ-Nglyco-Y_3 Y    
#> 10 20241224-LXJ-Nglyco-C_1 C    
#> 11 20241224-LXJ-Nglyco-C_2 C    
#> 12 20241224-LXJ-Nglyco-C_3 C
```

**📊 Expression Matrix - Your Data’s Heart and Soul:**

``` r
get_expr_mat(exp)[1:5, ]
#>                                       20241224-LXJ-Nglyco-H_1
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                  31054.12
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1                      NA
#> P04196-N344-Hex(5)HexNAc(4)                                NA
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               285613.66
#> P10909-N291-Hex(6)HexNAc(5)                       27588555.39
#>                                       20241224-LXJ-Nglyco-H_2
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                        NA
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1               136556.05
#> P04196-N344-Hex(5)HexNAc(4)                          15717.87
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               268250.71
#> P10909-N291-Hex(6)HexNAc(5)                       19527065.26
#>                                       20241224-LXJ-Nglyco-H_3
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                  457398.3
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1                      NA
#> P04196-N344-Hex(5)HexNAc(4)                          427312.0
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               2621248.6
#> P10909-N291-Hex(6)HexNAc(5)                        32930089.6
#>                                       20241224-LXJ-Nglyco-M_1
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                   7616346
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1                22675686
#> P04196-N344-Hex(5)HexNAc(4)                          10813133
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               993255511
#> P10909-N291-Hex(6)HexNAc(5)                          32500720
#>                                       20241224-LXJ-Nglyco-M_2
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                   7391049
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1                16675442
#> P04196-N344-Hex(5)HexNAc(4)                           9746325
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2              1099069766
#> P10909-N291-Hex(6)HexNAc(5)                                NA
#>                                       20241224-LXJ-Nglyco-M_3
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                   6267864
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1               114423292
#> P04196-N344-Hex(5)HexNAc(4)                          25348175
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2              1106268049
#> P10909-N291-Hex(6)HexNAc(5)                          26346060
#>                                       20241224-LXJ-Nglyco-Y_1
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                  23059718
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1               115717950
#> P04196-N344-Hex(5)HexNAc(4)                          33607210
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               547660361
#> P10909-N291-Hex(6)HexNAc(5)                          21780632
#>                                       20241224-LXJ-Nglyco-Y_2
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                  15010885
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1                90594397
#> P04196-N344-Hex(5)HexNAc(4)                          21284262
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               753702172
#> P10909-N291-Hex(6)HexNAc(5)                          19862189
#>                                       20241224-LXJ-Nglyco-Y_3
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                    740942
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1                55977605
#> P04196-N344-Hex(5)HexNAc(4)                          28608146
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               556784303
#> P10909-N291-Hex(6)HexNAc(5)                          12764805
#>                                       20241224-LXJ-Nglyco-C_1
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                        NA
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1               305840564
#> P04196-N344-Hex(5)HexNAc(4)                          32885077
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               669332806
#> P10909-N291-Hex(6)HexNAc(5)                          25946392
#>                                       20241224-LXJ-Nglyco-C_2
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                        NA
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1               428631806
#> P04196-N344-Hex(5)HexNAc(4)                          35418588
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2               696696106
#> P10909-N291-Hex(6)HexNAc(5)                          18860878
#>                                       20241224-LXJ-Nglyco-C_3
#> P08185-N176-Hex(5)HexNAc(4)NeuAc(2)                  10655.62
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-1             16064212.10
#> P04196-N344-Hex(5)HexNAc(4)                        6372648.51
#> P04196-N344-Hex(5)HexNAc(4)NeuAc(1)-2            112287653.59
#> P10909-N291-Hex(6)HexNAc(5)                        4316119.03
```

**The magic?** All the essential information for downstream analysis has
been carefully extracted, cleaned, and packaged into a beautiful
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object. No more data archaeology—just clean, analysis-ready data! ✨

## Peek Under the Hood: The Magic Explained 🔧✨

Curious about what just happened? Let’s lift the hood and see the
sophisticated machinery that
[`read_pglyco3_pglycoquant()`](https://glycoverse.github.io/glyread/reference/read_pglyco3_pglycoquant.md)
runs behind the scenes. Spoiler alert: it’s doing **a lot** of heavy
lifting so you don’t have to!

**The 8-Step Data Transformation Dance:**

1.  **Smart File Reading**: Reads your data with intelligent column type
    detection—no more “everything is a character” surprises!
2.  **Column Cleaning Magic**: Extracts clean UniProt accessions from
    messy “Proteins” columns and handles all those pesky formatting
    inconsistencies.
3.  **Protein Inference Intelligence**: Runs a sophisticated “parsimony”
    algorithm to resolve protein assignment ambiguities—because biology
    is complicated, but your data doesn’t have to be!
4.  **📊 PSM-to-Glycopeptide Aggregation**: Intelligently combines
    PSM-level quantification into meaningful glycopeptide-level
    measurements.
5.  **Sample Information Validation**: Cross-checks sample names and
    ensures everything matches up perfectly.
6.  **🎯 Data Extraction**: Carefully extracts variable information and
    expression matrix while maintaining data integrity.
7.  **🧬 Glycan Parsing Power**: Leverages the amazing
    [glyparse](https://github.com/glycoverse/glyparse) package to parse
    glycan compositions and structures into proper data types.
8.  **Final Assembly**: Packages everything into a beautiful,
    analysis-ready
    [`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
    object.

**The result?** What used to take hours of custom scripting now happens
in seconds, with bulletproof reliability!

## Meet the Whole glyread Family!

One function down, but there’s so much more! `glyread` is like a Swiss
Army knife for glycoproteomics data—it speaks the language of virtually
every major identification and quantification workflow out there. 🔧

**🌟 The Complete Toolkit:**

- **[`read_byonic_byologic()`](https://glycoverse.github.io/glyread/reference/read_byonic_byologic.md)**:
  The dynamic duo of Byonic identification + Byologic quantification
- **[`read_byonic_pglycoquant()`](https://glycoverse.github.io/glyread/reference/read_byonic_pglycoquant.md)**:
  Byonic’s precision meets pGlycoQuant’s power
- **[`read_msfragger()`](https://glycoverse.github.io/glyread/reference/read_msfragger.md)**:
  MSFragger-Glyco’s all-in-one identification and quantification magic
- **[`read_pglyco3()`](https://glycoverse.github.io/glyread/reference/read_pglyco3.md)**:
  Pure pGlyco3 with its built-in quantification capabilities
- **[`read_pglyco3_pglycoquant()`](https://glycoverse.github.io/glyread/reference/read_pglyco3_pglycoquant.md)**:
  The workflow we just explored—pGlyco3 + pGlycoQuant perfection
- **[`read_glyco_decipher()`](https://glycoverse.github.io/glyread/reference/read_glyco_decipher.md)**:
  Glyco-Decipher’s identification and label-free quantification
- **[`read_strucgp()`](https://glycoverse.github.io/glyread/reference/read_strucgp.md)**:
  StrucGP’s identification

**Coming Soon to a Package Near You:** We’re actively working on
supporting GPQuest, GlycanFinder and other exciting tools. The
glycoproteomics software landscape is evolving rapidly, and `glyread` is
evolving right alongside it! 🚀

**💡 Pro Tip:** No matter which workflow you choose, you can expect the
same consistent, clean data format on the other side. It’s like having a
universal remote for your glycoproteomics data!

## The Universal Data Language: Consistent Columns Across All Functions

Here’s the beautiful part: no matter which `glyread` function you use,
you’ll always get the same consistent, predictable data structure. It’s
like having a universal translator that ensures everyone speaks the same
language!

**🏷️ Your Variable Information Columns - The Standard Cast:**

- **`peptide`** 🧬: The peptide sequence (character) - your protein’s
  building blocks
- **`peptide_site`**: Glycosylation site on the peptide (integer) -
  where the magic happens
- **`protein`** 🏷️: UniProt accession (character) - your protein’s
  official ID card
- **`protein_site`** 🎯: Glycosylation site on the protein (integer) -
  the full-length coordinates
- **`gene`** 🧬: Gene symbol (character) - the genetic blueprint behind
  it all
- **`glycan_composition`**: A
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  object - what’s in your glycan
- **`glycan_structure`**: A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  object - how it’s all connected

**🎯 Special Note on Those Last Two:** The `glycan_composition` and
`glycan_structure` columns aren’t just plain text—they’re sophisticated
data types from the amazing
[glyrepr](https://github.com/glycoverse/glyrepr) package. Think of them
as smart objects that know how to do glycan math!

**💡 Pro Tip:** Getting familiar with `glyrepr` objects is absolutely
worth the investment. They unlock a whole world of glycan analysis
capabilities that would be impossible with plain text. Trust us on this
one!

## Your Glycoverse Adventure Awaits! 🚀✨

Congratulations! You’ve just mastered the art of taming glycoproteomics
data chaos with `glyread`. But this is just the beginning of your
journey through the `glycoverse` ecosystem!

**🎯 Your Next Destinations:**

- **[glyclean](https://github.com/glycoverse/glyclean)**: Your data
  preprocessing powerhouse! Think of it as Marie Kondo for
  glycoproteomics data—it’ll help you clean, filter, and organize your
  experiments until they spark joy.
- **[glymotif](https://github.com/glycoverse/glymotif)**: The motif
  hunting expedition! Discover hidden patterns and recurring structural
  themes in your glycan data. It’s like being a detective, but for sugar
  molecules!
- **📊 [glyexp](https://github.com/glycoverse/glyexp)**: Your
  experimental data command center! Master the
  [`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object and unlock the full power of synchronized data manipulation.

**The Best Part?** Since you’re now fluent in
[`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects thanks to `glyread`, all these packages will feel like natural
extensions of your workflow. It’s like learning a new language and
suddenly being able to read an entire library!

**Happy glycan hunting!** 🧬🎯

------------------------------------------------------------------------

*Remember: Great glycoproteomics analysis starts with great data import.
You’ve got the tools—now go make some discoveries! 💫*
