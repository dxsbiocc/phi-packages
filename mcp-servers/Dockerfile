FROM node:22-bookworm-slim

WORKDIR /app
COPY package.json package-lock.json ./
COPY mcp-servers ./mcp-servers
RUN npm ci && npm run build:mcp && npm prune --omit=dev

ENV HOST=0.0.0.0 \
    PORT=8080 \
    NODE_ENV=production \
    PHI_OMICS_CACHE_DIR=/var/cache/phi-omics-mcp
RUN mkdir -p /var/cache/phi-omics-mcp && chown -R node:node /var/cache/phi-omics-mcp /app
EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:'+(process.env.PORT||8080)+'/health').then(r=>process.exit(r.ok?0:1)).catch(()=>process.exit(1))"

USER node
CMD ["node", "mcp-servers/gateway/dist/cli.js"]
