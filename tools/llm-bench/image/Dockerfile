# Local-model benchmark image (#2095): the CI runner image (the game's build tools, EGL for scenarios) plus Aider, the
# coding agent the benchmark drives against the llama.cpp endpoint.
FROM raylibgame-runner:2.338.2
USER root
RUN apt-get update -q && apt-get install -y -q --no-install-recommends python3-venv python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/aider \
    && /opt/aider/bin/pip install -q --no-cache-dir aider-chat \
    && /opt/aider/bin/aider --version > /opt/aider/VERSION
ENV PATH=/opt/aider/bin:$PATH
USER runner
