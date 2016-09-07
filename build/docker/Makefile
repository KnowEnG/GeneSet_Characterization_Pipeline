CUR_DATE=`date +'%m_%d_%Y'`

.PHONY: build_docker_image

build_docker_image:
	docker build -t knowengdev/geneset_characterization:$(CUR_DATE) .

push_to_dockerhub:
	docker push knowengdev/geneset_characterization 
