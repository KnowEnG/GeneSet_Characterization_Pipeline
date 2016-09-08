CUR_DATE=`date +'%m_%d_%Y'`

build_docker_image:
	docker build -t knowengdev/geneset_characterization_pipeline:$(CUR_DATE) .

login_to_dockerhub:
	docker login --username=$(username) --email=$(email)

push_to_dockerhub:
	docker push knowengdev/geneset_characterization_pipeline
