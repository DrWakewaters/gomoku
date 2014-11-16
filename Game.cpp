#include "Game.h"

Game::Game(int maxSearchDepth, int boardHeight, int boardWidth, double maxSearchTime, bool blackIsHuman, bool whiteIsHuman, bool printAIInformation) : board(new Board(boardHeight, boardWidth)), blackPlayer(new Player(blackIsHuman, true)), whitePlayer(new Player(whiteIsHuman, false)), maxSearchDepth(maxSearchDepth), numberOfPliesMade(0), maxSearchTime(maxSearchTime), printAIInformation(printAIInformation) {}

Game::~Game() {
	delete this->board;
	delete this->blackPlayer;
	delete this->whitePlayer;
}

bool Game::makeBlackPly() {
	if(this->blackPlayer->isHuman) {
		return this->makeHumanPly(true);
	}
	return this->makeAIPly(true);
}

bool Game::makeWhitePly() {
	if(this->whitePlayer->isHuman) {
		return this->makeHumanPly(false);
	}
	return this->makeAIPly(false);
}

bool Game::makeHumanPly(bool isBlack) {
	Position position = {0, 0};
	bool empty = false;
	int y;
	int x;
	while(!empty) {
		std::cout << "Make a ply: ";
		std::cin >> y;
		std::cin >> x;
		position.y = static_cast<uint8_t>(y);
		position.x = static_cast<uint8_t>(x);
		if(!this->board->isOutsideBoard(position)) {
			empty = this->board->emptyAtPosition(position);
		}
	}
	if(isBlack) {
		this->board->setBlackAtPosition(position);
	} else {
		this->board->setWhiteAtPosition(position);		
	}
	this->board->updateInteresting(position);
	this->board->printBoard();
	this->board->clearHashTable();
	this->numberOfPliesMade++;
	this->board->updateSearchRectangle(position);
	return AI::gameWon(this->board, isBlack);
}

bool Game::makeAIPly(bool isBlack) {
	Position position;
	if(this->numberOfPliesMade > 0) {
		position = AI::computeScore(this->board, this->maxSearchDepth, this->maxSearchTime, isBlack, this->printAIInformation);
	} else {
		position = this->findStartingAIPly();
	}
	std::cout << "The ply chosen is (y, x) = (" << static_cast<int>(position.y) << ", " << static_cast<int>(position.x) << ")." << std::endl;
	if(isBlack) {
		this->board->setBlackAtPosition(position);
	} else {
		this->board->setWhiteAtPosition(position);
	}
	this->board->updateInteresting(position);	
	this->board->printBoard();
	this->board->clearHashTable();
	this->numberOfPliesMade++;
	this->board->updateSearchRectangle(position);
	return AI::gameWon(this->board, isBlack);
}

// If an AI player makes the first ply, let it put it somewhere in the middle of the board.
Position Game::findStartingAIPly() {
	return Position(9, 9);
	int8_t minY;
	int8_t maxY;
	int8_t minX;
	int8_t maxX;
	unsigned int seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	std::default_random_engine generator(seed);
	if(this->board->getBoardHeight() < 10) {
		minY = this->board->getBoardHeight()/2;
		maxY = this->board->getBoardHeight()/2;
	} else {
		minY = this->board->getBoardHeight()/2-3;
		maxY = this->board->getBoardHeight()/2+3;
	}
	if(this->board->getBoardWidth() < 10) {
		minX = this->board->getBoardWidth()/2;
		maxX = this->board->getBoardWidth()/2;
	} else {
		minX = this->board->getBoardWidth()/2-3;
		maxX = this->board->getBoardWidth()/2+3;
	}
	std::uniform_int_distribution<int8_t> yDistribution(minY, maxY);
	std::uniform_int_distribution<int8_t> xDistribution(minY, maxY);
	int8_t y = yDistribution(generator);
	int8_t x = xDistribution(generator);
	return Position(y, x);
}

bool Game::boardIsFull() {
	return this->board->boardIsFull();
}

int Game::getNumberOfPliesMade() {
	return this->numberOfPliesMade;
}

void Game::printBoard() {
	this->board->printBoard();
}
