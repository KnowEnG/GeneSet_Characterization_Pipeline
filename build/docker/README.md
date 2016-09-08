# Gene Set Characterization Pipeline
This Dockefile contains all the commands, in order, needed to build Gene Set Characterization pipeline. 

## Getting Started
Simply run the following command to build an image with latest code change. It uses the Dockerfile in current directory 
and generates a docker image with the tag indicating the date when it was created.
```makefile
    make build_docker_image
```
Then login to docker hub before you push to it. When prompted, enter your password and press enter.
```makefile
    make login_to_docker username=abc email=abc
```
Last, upload your image to docker hub!
```makefile
    make push_to_dockerhub
```


