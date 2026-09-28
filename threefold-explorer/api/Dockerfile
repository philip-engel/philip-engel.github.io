FROM sagemath/sagemath:10.9

USER root
WORKDIR /app
COPY engine/*.py /app/engine/
COPY engine/local_model_database.tar.gz /tmp/local_model_database.tar.gz
COPY service.py sage_worker.py /app/
RUN tar -xzf /tmp/local_model_database.tar.gz -C /app/engine \
    && rm /tmp/local_model_database.tar.gz \
    && chown -R sage:sage /app

USER sage
ENV PORT=8000 \
    PYTHONUNBUFFERED=1
EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=90s --retries=3 \
  CMD python3 -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/health', timeout=5)" || exit 1

CMD ["sage", "-python", "/app/service.py"]
