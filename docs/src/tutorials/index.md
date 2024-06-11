```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
  {
    href: "beginner/LoadingData", 
    src: "../vis.png",
    caption: "Loading Data with Pyehtim",
    desc: "How to load data using standard eht-imaging in Julia." 
  },
  {
    href: "beginner/GeometricModeling", 
    src: "../geom_model.png",
    caption: "Geometric Modeling of M87*",
    desc: "Modeling a black hole with simple geometric models" 
  }
];

const intermediate = [
  {
    href: "intermediate/ClosureImaging", 
    src: "../closure.png",
    caption: "Closure Imaging of M87*",
    desc: "Creating an image of a black hole using only closure information" 
  },
  {
    href: "intermediate/StokesIImaging", 
    src: "../stokesI.png",
    caption: "Simultaneous Imaging and Gain Modeling of M87*",
    desc: "Imaging a black hole with simultaneous gain modeling (selfcal) using complex visibilities" 
  },
  {
    href: "intermediate/PolarizedImaging", 
    src: "../telescopes.png",
    caption: "Full Stokes Imaging using RIME",
    desc: "Simultaneous instrument and polarized imaging of VLBI data." 
  }
];

const advanced = [
  {
    href: "advanced/HybridImaging", 
    src: "../hybrid.png",
    caption: "Hybrid ring modeling and residual imaging of M87*",
    desc: "How to combine everything to model the ring and create a residual image of M87*." 
  }
];


</script>



# Tutorials

This page contains a collection of tutorials that cover a range of topics from beginner to advanced. These demonstrate how to use Comrade in a variety of scenarios. While most of them consider the EHT, they should work more generally for any VLBI arrays.

## Beginner Tutorials

<Gallery :images="beginner" />

## Intermediate Tutorials

<Gallery :images="intermediate" />

## Advanced Tutorials

<Gallery :images="advanced" />


```
