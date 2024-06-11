```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
  {
    href: "beginner/LoadingData", 
    src: "../assets/vis.png",
    caption: "Loading Data with Pyehtim",
    desc: "How to load data using standard eht-imaging in Julia." 
  },
  {
    href: "beginner/GeometricModeling", 
    src: "../assets/geom_model.png",
    caption: "Geometric Modeling of M87*",
    desc: "Modeling a black hole with simple geometric models" 
  }
];

const intermediate = [
  {
    href: "intermediate/ClosureImaging", 
    src: "../assets/closure.png",
    caption: "Closure Imaging of M87*",
    desc: "Creating an image of a black hole using only closure information" 
  }
];


</script>



# Tutorials

## Beginner Tutorials

<Gallery :images="beginner" />

## Intermediate Tutorials

<Gallery :images="intermediate" />


```
