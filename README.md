# Emergence and widespread circulation of a recombinant SARS-CoV-2 lineage in North America
Investigation of the potential role of recombination on the origin and evolution of SARS-CoV-2 lineages B.1.627, B.1.628, B.1.631 and B.1.634.


Bernardo Gutierrez<sup>1,2,3</sup>, Hugo Gildardo Castelán Sánchez<sup>1,2,4</sup>, Darlan Cândido da Silva<sup>1,5</sup>, Ben Jackson<sup>6</sup>, Shay Fleishon<sup>7</sup>, Christopher Ruis<sup>8,9</sup>, Luis José Delaye Arredondo<sup>2,10</sup>, Andrew Rambaut<sup>6</sup>, Oliver G. Pybus<sup>1,11*</sup>, Marina Escalera-Zamudio<sup>1,2*</sup>


<sup><sup>1</sup>Department of Zoology, University of Oxford, Oxford, UK.
<sup>3</sup>Consorcio Mexicano de Vigilancia Genómica (CoViGen-Mex).
<sup>3</sup>Colegio de Ciencias Biológicas y Ambientales, Universidad San Francisco de Quito, Quito, Ecuador.
<sup>3</sup>Instituto de Microbiología, Universidad San Francisco de Quito, Quito, Ecuador.
<sup>4</sup>CONACYT-UAEM-Centro de Investigación en Dinámica Celular, Mexico.
<sup>5</sup>Instituto de Medicina Tropical, Faculdade de Medicina da Universidade de São Paulo, Brazil.
<sup>6</sup>Institute of Evolutionary Biology, University of Edinburgh, UK.
<sup>7</sup>Central Virology Laboratory, Public Health Services, Ministry of Health, Sheba Medical Center, Tel-Hashomer 52621, Israel.
<sup>8</sup>Molecular Immunity Unit, Department of Medicine, University of Cambridge, Cambridge, UK.
<sup>9</sup>Department of Veterinary Medicine, University of Cambridge, Cambridge, UK.
<sup>10</sup>Departamento de Ingeniería Genética, Unidad Irapuato, CINVESTAV, Mexico.
<sup>11</sup>Department of Pathobiology, Royal Veterinary College, London UK.</sup>


## Repository usage and structure

The structure of this repository is shown below. All R scripts used to generate plots or analyses are located in the root directory. Further details for the contents of specific directories is further explained below.

```
XB_lineage_investigation/
├── Data
│   ├── OWID_data
│   └── genetic_data
│       ├── GISAID_20210831_alignments
│       ├── Reduced_alignments
│       ├── Snipit_alignments
│       └── metadata
├── XB_phylogenetics
│   ├── IQ_Tree
│   │   ├── GISAID_20210831_trees
│   │   └── Reduced_trees
│   └── UShER
├── XB_recombination_analysis
│   ├── 3SEQ
│   ├── RDP
│   └── GARD
├── Plots
├── Acknowledgements
└── README.md
```

## Input data

The [`Data`](Data/) directory contains specific alignments used for the phylogenetic and recombination analyses. All SARS-CoV-2 genomic data was retrieved from GISAID (www.gisaid.org). We are grateful to all the laboratories and institutions worldwide involved in the generation of virus genome data shared on GISAID, and a full list of the laboratories who produced the sequences used in this study is available in the [`Acknowledgements`](Acknowledgements/) directory. Additional epidemiological data was retrieved from the Our Wrld in Data repository (https://ourworldindata.org/coronavirus).
