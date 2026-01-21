# MetalinksDB Database

Pre-compiled metabolite-protein interaction database from MetalinksDB

## Usage

``` r
metalinksdb
```

## Format

A list containing the following components:

- edges:

  data.frame with 41,894 metabolite-protein interactions containing:

  - hmdb: HMDB metabolite identifier

  - uniprot: UniProt protein identifier

  - source: Data source

  - combined_score: Interaction confidence score (0-1000)

  - mor: Mode of regulation (1=producing, -1=degrading, 0=binding)

  - type: Interaction type ("lr"=ligand-receptor, "pd"=produce-degrade)

  - transport_direction: For transporters ("in" or "out")

- metabolites:

  data.frame with 1,128 metabolites containing:

  - hmdb: HMDB identifier

  - metabolite: Metabolite name

  - pubchem: PubChem identifier

  - metabolite_subclass: Chemical classification

- proteins:

  data.frame with 4,374 proteins containing:

  - uniprot: UniProt identifier

  - gene_symbol: Gene symbol

  - protein_type: Protein classification (enzyme, gpcr, transporter,
    etc.)

- pathway:

  data.frame with 157,741 metabolite-pathway associations

- cell_location:

  data.frame with 2,816 subcellular location annotations

- tissue_location:

  data.frame with 2,410 tissue location annotations

- disease:

  data.frame with 3,216 disease associations

## Source

MetalinksDB (https://metalinks.org/)

## Details

The database enables two types of metabolite-mediated communication
inference:

\*\*Ligand-Receptor (lr) type\*\*: Direct metabolite-receptor binding
interactions, primarily involving GPCRs, nuclear hormone receptors, and
ion channels.

\*\*Produce-Degrade (pd) type\*\*: Enzyme-mediated metabolite production
and consumption, enabling inference of metabolite availability from
enzyme expression.

## References

Schafer, S., et al. (2023). MetalinksDB: a knowledgebase of
metabolite-centric signaling. Nature Communications.

## Examples

``` r
# \donttest{
# Access the database
db <- scMetaLink:::.load_metalinksdb()

# View available metabolites
head(db$metabolites)
#> # A tibble: 6 × 4
#>   hmdb        metabolite        pubchem metabolite_subclass                     
#>   <chr>       <chr>             <chr>   <chr>                                   
#> 1 HMDB0001448 Sulfate           1117    Non-metal sulfates                      
#> 2 HMDB0000208 Oxoglutaric acid  51      Gamma-keto acids and derivatives        
#> 3 HMDB0000464 Calcium           271     NA                                      
#> 4 HMDB0000131 Glycerol          753     Carbohydrates and carbohydrate conjugat…
#> 5 HMDB0000429 17alpha-Estradiol 68570   Estrane steroids                        
#> 6 HMDB0000990 Acetaldehyde      177     Carbonyl compounds                      

# Check interaction types
table(db$edges$type)
#> 
#>    lr    pd 
#> 11869 30025 
# }
```
