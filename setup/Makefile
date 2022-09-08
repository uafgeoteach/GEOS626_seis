all: build run

build:
	docker build -f dockerfile -t geos626_seis:latest .

build-no-cache:
	docker build --pull --no-cache -f dockerfile -t geos626_seis:latest .

run:
	bash run.sh 2>&1 | tee log

clean:
	rm -rf ../.cache ../.ipython .../jupyter .../local \
		../.bash_history ../.ipynb_checkpoints ../__pycache__
