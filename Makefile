test:
	bash etc/test.sh

build_docker:
	bash etc/build_docker.sh

push_docker:
	bash etc/build_docker.sh push_yes

.PHONY: test build_docker
