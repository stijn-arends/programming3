<div id="top"></div>

<!-- [![Contributors][contributors-shield]][contributors-url] -->
<!-- [![MIT License][license-shield]][license-url] -->

<div align="center">
<h1 align="center">Assignment 6</h1>
<h2 align="center">Investigation of the scientific literature</h2>
  <!-- <a href="https://www.molgenis.org/">
    <img src="figures/PubMed-Logo.png" alt="Logo" width="100" height="100">
  </a>
  <a>
    <img src="figures/Apache_Spark_logo.png" alt="spark" width="100" height="100" align='right'>
  </a>
   <a>
    <img src="figures/networkx_logo.svg" alt="spark" width="100" height="100" align='left'>
  </a> -->
</div>
<div class="row">
  <div class="column" style='float:left;width:33.33%;padding:5px'>
    <img src="figures/networkx_logo.svg" alt="pubmed" style="width:100%">
  </div>
  <div class="column" style='float:left;width:33.33%;padding:5px'>
    <img src="figures/PubMed-Logo.png" alt="Forest" style="width:100%">
  </div>
  <div class="column" style='float:left;width:33.33%;padding:5px'>
    <img src="figures/Apache_Spark_logo.png" alt="Mountains" style="width:100%">
  </div>
</div>

# Assignment 6


## Getting Started

Server:
```bash
python src/main.py -d /data/dataprocessing/NCBI/PubMed/ -p 4235 --host assemblix2019 -s -o /commons/dsls/dsph/2022/test/
```

Client: 
```bash
python src/main.py -p 4235 --host assemblix2019 -c -n 2
```

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/stijn-arends/programming3.svg?style=for-the-badge
[contributors-url]: https://github.com/stijn-arends/programming3/graphs/contributors
[license-shield]: https://img.shields.io/github/license/stijn-arends/programming3.svg?style=for-the-badge
[license-url]: https://github.com/stijn-arends/programming3/blob/master/LICENSE.md
