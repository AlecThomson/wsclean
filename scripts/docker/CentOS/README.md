# Introduction

Dockerfile to build a WSClean docker image (master branch) based on CentOS.


# How to build the docker image

Just type the following commands

```
$ make  # for build/release version
#or
$ make build # just for build
#or
$ make release # this will generate a version base on master file
#or
$ make push # this way will make other people can pull the docker image
```