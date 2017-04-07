from celery import Celery

app = Celery('task', broker='amqp://localhost//')

@app.task
def reverse(name):
    return name[::-1]