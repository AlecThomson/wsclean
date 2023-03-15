OS=centos8-wsclean

all: build release

VERSION='master'

build:
	@docker build -f Centos8 --tag=aroffringa/$(OS):latest .

release: build
	@docker build -f Centos8 --tag=aroffringa/$(OS):master$(VERSION) .

push:
	@docker push aroffringa/$(OS):$(VERSION)
	@docker push aroffringa/$(OS):latest