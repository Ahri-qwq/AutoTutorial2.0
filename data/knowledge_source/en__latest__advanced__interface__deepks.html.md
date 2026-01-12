# https://abacus.deepmodeling.com/en/latest/advanced/interface/deepks.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/interface/deepks.html

# DeePKS[#](#deepks)

[DeePKS](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00872) is a machine-learning (ML) aided density funcitonal model that fits the energy difference between highly accurate but computationally demanding method and effcient but less accurate method via neural-network. Common high-precision methods include hybrid functionals or CCSD-T, while common low-precision methods are LDA/GGA.

[J. Phys. Chem. A 126.49 (2022): 9154-9164](https://pubs.acs.org/doi/abs/10.1021/acs.jpca.2c05000) for the application of DeePKS in periodic systems.

[website](https://deepks-kit.readthedocs.io/en/latest/index.html). The DeePKS-related keywords in `INPUT`

file can be found [here](http://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#deepks). An [example](https://github.com/deepmodeling/deepks-kit/tree/abacus/examples/water_single_lda2pbe_abacus) for training DeePKS model with ABACUS is also provided. For practical applications, users can refer to a series of [Notebooks](https://bohrium.dp.tech/collections/1921409690). These Notebooks provide detailed instructions on how to train and use the DeePKS model using perovskite as an example. Currently, these tutorials are available in Chinese, but we plan to release corresponding English versions in the near future.

Note: DeePKS calculations can only be performed by the LCAO basis.
