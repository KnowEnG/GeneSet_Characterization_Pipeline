CUR_DATE=`date +'%m_%d_%Y'`

all: build_docker_image login_to_dockerhub push_to_dockerhub

build_docker_image:
	docker build --no-cache -t knowengdev/geneset_characterization_pipeline:$(CUR_DATE) .

login_to_dockerhub:
	docker login --username=$(username) --password=$(password)

push_to_dockerhub:
	docker push knowengdev/geneset_characterization_pipeline:$(CUR_DATE)
