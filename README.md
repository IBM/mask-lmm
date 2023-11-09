<!-- ABOUT THE PROJECT -->

## :dna: MaSK-LMM: Matrix-Sketching Linear Mixed Model

[![MaSK-LMM Diagram][masklmm]](#)

<!-- <a>![badge-alt-text](images/comical.jpg)</a> -->

We leverage _matrix sketching_ to develop a fast and efficient LMM method called **Ma**trix-**Sk**etching **LMM** (**MaSk-LMM**) by sketching the genotype matrix to reduce its dimensions and speed up computations.

### Built With

<!-- This section should list any major frameworks/libraries used to bootstrap your project. Leave any add-ons/plugins for the acknowledgements section. Here are a few examples. -->

[![Anaconda][Anaconda.com]][Anaconda-url]
[![Python][Python.com]][Python-url]
[![VSCode][VSCode.com]][VSCode-url]
[![Github][Github.com]][Github-url]

<!-- [![Jupyter][Jupyter.com]][Jupyter-url] -->

## Getting Started
<!--
### Prerequisites

1. Create the environment from the `comical_env.yml` file:

```

conda env create -f comical_env.yml

```

* Note: if you receive the error `bash: conda: command not found...`, you need to install Anaconda to your development environment (see "Additional resources" below)

2. Activate the new environment:

```

conda activate comical-env

```

3. Verify that the new environment was installed correctly:

```

conda env list

```

* Additional resources:

   * [Connect to computing cluster](http://ccc.pok.ibm.com:1313/gettingstarted/newusers/connecting/)

   * [Set up / install Anaconda on remote linux server](https://kengchichang.com/post/conda-linux/)

   * [Set up remote development environment using VSCode](https://code.visualstudio.com/docs/remote/ssh)

<a name="running_comical"></a>
-->

### Installing MaSK-LMM

```
pip install masklmm
```

### Running MaSK-LMM


```python
# Data

data_path = "../sample-data/bn"
bed_fn = data_path+".bed"
pheno_fn = data_path+".phen_w_header"
cov_fn = data_path+".cov_w_header"
pruned_bed_fn = "../sample-data/bn.bed"

# Parameters

maxiters = 10
sample_sketch_size = 0.5
marker_sketch_size = 0.5
block_size = 10000

# MaSkLMM pacakge

from masklmm import MaSkLMM

# Running MaSkLMM

newton = MaSkLMM.compute(bed_fn,
                        pruned_bed_fn,
                        pheno_fn,
                        cov_fn,
                        sample_sketch_size = sample_sketch_size,
                        marker_sketch_size = marker_sketch_size,
                        maxiters = maxiters,
                        block_size = block_size)
```


<!--
### Help

```

python wrapper.py --help

```
-->
## Authors

Contributors and contact info:

* Myson Burch (myson dot burch at ibm dot com)

* Aritra Bose (a dot bose at ibm dot com)


## Version History

<!-- * 0.2

    * Various bug fixes and optimizations

    * See [commit change]() or See [release history]() -->

* 0.1.0

    * Initial Release

<!-- ## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details -->

<!-- MARKDOWN LINKS & IMAGES -->

<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge

[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors

[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge

[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members

[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge

[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers

[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge

[issues-url]: https://github.com/othneildrew/Best-README-Template/issues

[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge

[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt

[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555

[linkedin-url]: https://linkedin.com/in/othneildrew

[masklmm]: images/masklmm_overview.png

[notebook]: images/screenshot.png

[VSCode.com]: https://img.shields.io/badge/Visual_Studio_Code-033b66?style=for-the-badge&logo=visual%20studio%20code&logoColor=white

[VSCode-url]: https://code.visualstudio.com

[Python.com]: https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54

[Python-url]: https://www.python.org

[Jupyter.com]: https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white 

[Jupyter-url]: https://jupyter.org

[Github.com]: https://img.shields.io/badge/github-%23006567.svg?style=for-the-badge&logo=github&logoColor=white 

[Github-url]: https://github.com

[Anaconda.com]: https://img.shields.io/badge/Anaconda-%2344A833.svg?style=for-the-badge&logo=anaconda&logoColor=white

[Anaconda-url]: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment


