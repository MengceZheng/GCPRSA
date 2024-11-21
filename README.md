# Generalized Cryptanalysis of Cubic Pell RSA

## Introduction

This is a Python implementation of lattice-based attack proposed in **Generalized Cryptanalysis of Cubic Pell RSA**[^GCPRSA].

## Requirements

- [**SageMath**](https://www.sagemath.org/) 9.5 with Python 3.10

You can check your SageMath Python version using the following command:

```commandline
$ sage -python --version
Python 3.10.12
```

Note: If your SageMath Python version is older than 3.9.0, some features in given scripts might not work.

## Usage

The standard way to run the attack with given parameters $N$, $e$, $\gamma$, $m$ and $t$ requires passing them as arguments step by step in command line. For instance, to run the attack with given parameters used in our paper, please first run `sage -python attack.py`:

```commandline
GCPRSA$ sage -python attack.py
Input given parameters of GCPRSA attack instance as follows
Input N: 550366209463983254224851898151920438687572141757121552287270257270437967965957081683577937037276073506051924501113396260170171
Input e: 1057800388414613269699393034579590821261008821244344311613132208333483522085457259293084165274518494991109201666203006752031456045032171612863063434022522609550692892561154763861484988711873034869148741612190479043963664788377209
Input g: 0.5
Input m: 4
Input t: 1
Found primes: p = 967502495361032247552444598347042412041475154993790090306919213 and q = 568852496094709460190033647617209776684578742299012948180863367
The attack costs 0.846 seconds...
```

## Notes

All the details of the numerical attack experiments are recorded in the `attack.log` file.

[^GCPRSA]: Kang H., Zheng M., "Generalized Cryptanalysis of Cubic Pell RSA"
