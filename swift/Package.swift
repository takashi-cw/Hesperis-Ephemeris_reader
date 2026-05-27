// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "StellaBspReader",
    platforms: [
        .macOS(.v14),
    ],
    products: [
        .library(
            name: "StellaBspReader",
            targets: ["StellaBspReader"]
        ),
    ],
    targets: [
        .target(
            name: "StellaBspReader",
            path: "Sources/StellaBspReader"
        ),
        .testTarget(
            name: "StellaBspReaderTests",
            dependencies: ["StellaBspReader"],
            path: "Tests/StellaBspReaderTests"
        ),
    ]
)
