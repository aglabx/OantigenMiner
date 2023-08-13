
<a href=""><img src="img/woof_logo.png" align="right" width="200" ></a>

# WOOF 


### Wonderful O-antigens Operons Finder

<br />

> **Authors:** <br />
Oksana Kotovskaia <br />
Ekaterina Marenina <br />
Igor Ostanin <br />
Nadezhda Pavlova <br />
Nikita Vaulin <br />
*Independent researchers from Glasgow* <br /><br />
**Supervisor:** <br />
Polina Kuchur, [*chesnokova@scamt-itmo.ru*](mailto:chesnokova@scamt-itmo.ru), <br /> *Applied genomics laboratory, ITMO University, Saint-Petersburg, Russia*

## Overview

O-antigens are essential components of lipopolysaccharides exposed on the surface of bacterial cells, and their diversity has increased significantly over time. The identification of O-antigen operons is crucial for comparative genomics, bacterial identification, epidemiological studies, and vaccine development tasks. Here we introduce **WOOF**, an automatic pipeline for searching and visualizing O-antigen operons.

 
## Installation

You can get **WOOF** *via* GitHub

```
git clone git@github.com:aglabx/OantigenMiner.git
```

Then move to the repository directory:

```
cd OantigenMiner
```

## Usage

### Example usage

Write path to fasta file:
```
data = Path({folder_name})
GENOME = "{file_prefix}"
email = 'myemail@yarevan-hackaton.io'
```
Run WOOF:
```
snakemake -c {threads} --use-conda
```


### Arguments


## Uninstallation


## Troubleshooting

If you encounter any issues or have questions, you can open an issue on the GitHub page or directly contact the authors for support.

## Citation

If you use these tool, please cite as:
- Kotovskaia O., Ostanin I., Marenina E., Pavlova N., Vaulin N. , Kuchur P., Komissarov A. WOOF: Wonderful O-anrtigens Operons Finder. Природа, 2023, 220, 1--128
```bibtex
@article{phandorin2023,
  title={WOOF: Wonderful O-anrtigens Operons Finder},
  author={Kotovskaia O., Ostanin I., Marenina E., Pavlova N., Vaulin N. , Kuchur P., Komissarov A.},
  journal={Природа},
  year={2023},
  volume={220},
  pages={1--128}
}
```

