# Building The GeneSet Characterization Pipeline Docker Image
The Dockefile in this directory contains all the commands, in order, needed to build the geneset_characterization_pipeline docker image.

* run the "make" command to build the samples_clustering_pipeline docker image (output: docker image called "geneset_characterization_pipeline" and a tag with today's date and time):
```
    make build_docker_image
```

* Login to docker hub. When prompted, enter your password and press enter:
```
    make login_to_docker username=your docker login here email=your email as in your Docker profile here
```

* Upload your image to docker hub:
```
    make push_to_dockerhub
```


