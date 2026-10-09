import fs from "node:fs";
import type { IncomingMessage, ServerResponse } from "node:http";
import path from "node:path";
import type { Plugin } from "vite";
import { defineConfig } from "vitest/config";

/**
 * Serves finished runs and slices from `../out` (or $WAVESIM_RUNS) under `/runs/<name>/...`, and
 * a listing of both at `/runs/index.json`. Runs are large and git-ignored, so they are never
 * copied into the build.
 */
function runs(dir: string): Plugin {
  const types: Record<string, string> = {
    ".json": "application/json",
    ".f32": "application/octet-stream",
  };

  const handle = (req: IncomingMessage, res: ServerResponse, next: () => void) => {
    const url = new URL(req.url ?? "/", "http://localhost");
    if (!url.pathname.startsWith("/runs/")) return next();
    const rel = decodeURIComponent(url.pathname.slice("/runs/".length));

    if (rel === "index.json") {
      const holding = (file: string): string[] =>
        fs.existsSync(dir)
          ? fs
              .readdirSync(dir, { withFileTypes: true })
              .filter((e) => e.isDirectory() && fs.existsSync(path.join(dir, e.name, file)))
              .map((e) => e.name)
              .sort()
          : [];
      res.setHeader("Content-Type", "application/json");
      res.end(JSON.stringify({ runs: holding("run.json"), slices: holding("slice.json") }));
      return;
    }

    const file = path.resolve(dir, rel);
    if (!file.startsWith(path.resolve(dir) + path.sep)) {
      res.statusCode = 403;
      res.end("forbidden");
      return;
    }
    if (!fs.existsSync(file) || !fs.statSync(file).isFile()) {
      res.statusCode = 404;
      res.end("not found");
      return;
    }
    res.setHeader("Content-Type", types[path.extname(file)] ?? "application/octet-stream");
    res.setHeader("Content-Length", fs.statSync(file).size);
    res.setHeader("Cache-Control", "no-cache");
    fs.createReadStream(file).pipe(res);
  };

  return {
    name: "wavesim-runs",
    configureServer: (server) => void server.middlewares.use(handle),
    configurePreviewServer: (server) => void server.middlewares.use(handle),
  };
}

export default defineConfig({
  plugins: [runs(path.resolve(process.env.WAVESIM_RUNS ?? path.join(import.meta.dirname, "../out")))],
  build: {
    target: "es2022",
    chunkSizeWarningLimit: 800,
    rollupOptions: {
      input: {
        main: path.join(import.meta.dirname, "index.html"),
        slice: path.join(import.meta.dirname, "slice.html"),
      },
    },
  },
  // Two test workers at most: the suite is small, and this keeps a laptop quiet.
  test: { include: ["tests/**/*.test.ts"], maxWorkers: 2 },
});
