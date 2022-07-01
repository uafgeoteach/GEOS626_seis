all: build run

build:
	docker build -f setup/dockerfile -t geos626_seis:latest .

run:
	bash run_geos626_container.sh 2>&1 | tee log

clean:
	rm -rf home/.cache home/.ipython home/.jupyter home/.local home/.bash_history