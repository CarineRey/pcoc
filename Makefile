test:
	bash etc/test.sh
test_det:
	bash etc/test_det.sh

build_docker:
	bash etc/build_docker.sh

push_docker:
	bash etc/build_docker.sh push_yes

.PHONY: test build_docker test_det
