# Gene Set Characterization Pipeline
This Dockefile contains the necessity of building the environment for running Gene Set Characterization pipeline. 

## Getting Started
Simply run the following command to build an images with latest code chage. It uses the Dockerfile in current directory 
and generates a docker image with the tag indicating the date of current day.
```
    make build_docker_image
```
Then login to docker hub before you push to it. When prompted, enter your password and press enter.
```
    make login_to_docker username=abc email=abc
```
Last, upload your image to docker hub!
```
    make push_to_dockerhub
```


