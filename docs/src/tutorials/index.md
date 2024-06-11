```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
  {
    href: "beginner/LoadingData", 
    caption: "Loading Data with Pyehtim", 
    desc: "How to load data using standard eht-imaging in Julia." 
  },
  {
    href: "beginner/GeometricModeling.md",
    caption: "Geometric Modeling of M87",
    desc: "Fitting simple geometric models to the 2017 M87 data."
  }
];

</script>



# Tutorials

## Beginner Tutorials

<Gallery :images="beginner" />


```
