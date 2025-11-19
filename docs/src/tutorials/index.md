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
    src: "../telescopes.jpg",
    caption: "Full Stokes Imaging using RIME",
    desc: "Simultaneous instrument and polarized imaging of VLBI data." 
  }
];

const advanced = [
  {
    href: "advanced/Hibi", 
    src: "../telescopes.jpg",
    caption: "Hierarchical Bayesian Imaging (HIBI) of M87*",
    desc: "How to incorporate physical assumptions in to the modeo." 
  },
  {
    href: "advanced/FitPS", 
    src: "../telescopes.jpg",
    caption: "Fitting the Power Spectrum of an AGN",
    desc: "Fitting a power spectrum model to VLBI data of an AGN." 
  },

];


</script>



# Tutorials

This page contains a collection of tutorials that cover a range of topics from beginner to advanced. These demonstrate how to use Comrade in a variety of scenarios. While most of them consider the EHT, they should work more generally for any VLBI arrays.

!!! warning
    All plots in this tutorial are shown using `DisplayAs` to 
    prevent the webpage from being too big and slow. You do not
    need to include `|> DisplayAs.PNG |> DisplayAs.Text` in your
    code when running locally.

## Beginner Tutorials

<Gallery :images="beginner" />

## Intermediate Tutorials

<Gallery :images="intermediate" />

## Advanced Tutorials

<Gallery :images="advanced" />


```
